#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    gel the clusters lists 

"""

import pandas as pd
import pickle
import collections
import scipy.stats as spicy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator




 ###################################### functions  ######################################


def make_enriched_df(enriched_list, name_df):
    '''
    create a new df with enriched genes (E or C) and pres/abs on all strains
    '''
    data = []
    for gene in enriched_list:
        data.append(df.loc[df['Gene']==gene])
    newdf = pd.DataFrame(columns=df.columns)
    newdf.name = name_df
    for line in data : 
        newdf = newdf.append(line)
    newdf = newdf.reset_index(drop=True)
    
    return(newdf)



 ###################################### main  ######################################


    ########## generate df of enriched genes with the position of each one in each strain ##########

gpa = pd.read_csv('../gene_infos/gpa.csv', header = 0,dtype=str) #gene_presence_absence.csv from roary

df = gpa.drop(['Pref_name','Func_COG','Description','Annotation','No_isolates'],axis=1)

dfC = pd.read_pickle('../pickle_stuff/enriched_C.pickle')
enriched_C = dfC['Gene'].tolist()
df_C = make_enriched_df(enriched_C, 'df_C')

dfE = pd.read_pickle('../pickle_stuff/enriched_E.pickle')
enriched_E = dfE['Gene'].tolist()
df_E = make_enriched_df(enriched_E, 'df_E')



    ########## generate clusters of enriched genes in C (group 1) ##########

C_strains = []                                                #list of tuples [({clust},C/E)]

for nb_col,colName in enumerate(df_C):
    
    if colName == 'Gene':
        continue

    cluster = set()                                             #reset par souche + à la fin de chaque bloc trouvé
    
    tmp_df = df_C.filter(['Gene',colName])
    tmp_df.dropna(axis=0,inplace=True)
    tmp_df.sort_values(colName,inplace=True)
    tmp_df.reset_index(inplace=True,drop=True) 
    
    for index,line in tmp_df.iterrows() :
        if pd.isnull(tmp_df[colName][index]) :                  #skip gène si pas présent dans la souche testée
            continue

        pos1 = int(tmp_df[colName][index].split('_')[1])        #prendre pos 1
        if index+1 in range(len(tmp_df)):                       #get le dernier gène 
            pos2 = int(tmp_df[colName][index+1].split('_')[1])  #prendre pos2

            if pos2 - pos1 <= 3 :                               #si mes deux pos (gènes) ont max 2 (gènes) d'écarts :
                cluster.add(tmp_df['Gene'][index])
                cluster.add(tmp_df['Gene'][index+1])            #add les 2 gènes au bloc en cours (set -> pas de doublons)

            else :                                              #si on est à la fin d'un cluster 
                if len(cluster) > 3 :                           #au moins 4 gènes     
                    tup=cluster,colName[0]                      #créer tuple set(cluster),E/C
                    C_strains.append(tup)                       #on l'ajoute à la liste de cluster asso à E/C selon colName
                cluster = set()                                 #et on reset le cluster en cours qu'il ait été add ou pas
                    
    if len(cluster) > 3 :                                       #si mon dernier cluster > 3 gènes
        tup=cluster,colName[0]                        
        C_strains.append(tup)                                   #on l'ajoute aussi (avec colName[0])




    #sort and deduplicate

sorted_C_strains = sorted(C_strains,key=lambda x : len(x[0]),reverse=True)

#je veux [ [{cluster} , nbE , nbC] ] -> list of list :)

deduplicated_C_strains = []

for i in range(len(sorted_C_strains)):
    block,strain = sorted_C_strains[i]   
    index = [x for x,y in enumerate(deduplicated_C_strains) if y[0]==block]
    
    if len(index)==0 :
        if strain == 'E':
            deduplicated_C_strains.append([block,1,0])
        elif strain == 'C': 
            deduplicated_C_strains.append([block,0,1])
    
    elif len(index)==1:
        index = index[0]
        group = deduplicated_C_strains[index]
        
        if strain == 'E':
            group[1]=group[1]+1

        elif strain == 'C':
            group[2]=group[2]+1



#   sort par occ du dernier chiffre == C, puis par len si same occ
sorted_deduplicated_C_strains = sorted(deduplicated_C_strains,key=lambda x:(x[-1],len(x[0])),reverse=True) 


C_to_group = []

for block,E,C in sorted_deduplicated_C_strains :
    if C >= 2:
        C_to_group.append((block,E,C))



    #group the clusters and determine if is subset or superset of an already existing cluster, keep the most occuring


grouped_C_strains = []

for i in range(len(C_to_group)):
    
    block,E,C = C_to_group[i]
    
    common = [element for element,E,C in grouped_C_strains if not block.isdisjoint(element)]
    
    if len(common)==0 :
        grouped_C_strains.append((block,E,C))
#        print('block',i+1,'( len',len(block),', occ E/C',E,C,') of C_to_group added to grouped list as block',grouped_C_strains.index((block,E,C))+1)

        
    elif len(common)==1:
        issubset = [element for element,E,C in grouped_C_strains if block.issubset(element)]
        issuperset = [element for element,E,C in grouped_C_strains if block.issuperset(element)]
        
        if len(issubset)==1:
            clust = issubset[0]
            clust_to_search = [tup for tup in grouped_C_strains if tup[0]==clust]
            block2,E2,C2 = clust_to_search[0]
#            print('block',i+1,'( len',len(block),', occ E/C',E,C,') of C_to_group is sub set of block', 
#                   grouped_C_strains.index((block2,E2,C2))+1,'( len',len(block2),', occ',E2,C2,') in grouped list')
       
#        elif len(issubset)>1:
#            print('block is subset of',len(issubset),'other groups in grouped list:\n',issubset)
        
           
        if len(issuperset)==1:
            clust = issuperset[0]
            clust_to_search = [tup for tup in grouped_C_strains if tup[0]==clust]
            block2,E2,C2 = clust_to_search[0]
#            print('block',i+1,'( len',len(block),', occ E/C',E,C,') of C_to_group is super set of block', 
#                 grouped_C_strains.index((block2,E2,C2))+1,'( len',len(block2),', occ E/C',E2,C2,') in grouped list')
        
        #elif len(issuperset)>1:
        #    print('block is superset of',len(issuperset),'other groups in grouped list:\n',issuperset)
                
    #else : 
    #    print('genes in block',i+1,'( len',len(block),'occ E/C',E,C,') commons with',len(common),'other group.s :\n',common)



        #write the results for group one : one excel with all, one excel with one cluster per sheet
row_list = []

for group,E,C in grouped_C_strains :
    for gene in group : 
        
        dico = {}
        dico['nb_clust'] = grouped_C_strains.index((group,E,C))+1
        dico['nb_gene']=len(group)
        dico['occ_E']=E
        dico['occ_C']=C
        dico['gene']=gene
        dico['pref_name']=gpa.loc[int(gpa[gpa['Gene']==gene].index.values),'Pref_name']
        dico['func']=gpa.loc[int(gpa[gpa['Gene']==gene].index.values),'Func_COG']
        dico['description']=gpa.loc[int(gpa[gpa['Gene']==gene].index.values),'Description']
        dico['annot']=gpa.loc[int(gpa[gpa['Gene']==gene].index.values),'Annotation']
        
        row_list.append(dico)
        
df_clust_C = pd.DataFrame(row_list)
df_clust_C.to_csv('../new_clust_C.csv',index=False)
writer = pd.ExcelWriter('new_clusters_C.xlsx')

for group,data in df_clust_C.groupby('nb_clust'):
    data.to_excel(writer,sheet_name = str(group),index=False)
writer.save()




    ########## generate clusters of enriched genes in E (group 2) ##########

E_strains = []                                                #list of tuples [({clust},C/E)]

for nb_col,colName in enumerate(df_E):
    
    if colName == 'Gene':
        continue

    cluster = set()                                             #reset par souche + à la fin de chaque bloc trouvé
    
    tmp_df = df_E.filter(['Gene',colName])
    tmp_df.dropna(axis=0,inplace=True)
    tmp_df.sort_values(colName,inplace=True)
    tmp_df.reset_index(inplace=True,drop=True) 
    
    for index,line in tmp_df.iterrows() :
        if pd.isnull(tmp_df[colName][index]) :                  #skip gène si pas présent dans la souche testée
            print('ya un pb :)')
            continue

        pos1 = int(tmp_df[colName][index].split('_')[1])        #prendre pos 1
        if index+1 in range(len(tmp_df)):                       #get le dernier gène 
            pos2 = int(tmp_df[colName][index+1].split('_')[1])  #prendre pos2

            if pos2 - pos1 <= 3 :                               #si mes deux pos (gènes) ont max 2 (gènes) d'écarts :
                cluster.add(tmp_df['Gene'][index])
                cluster.add(tmp_df['Gene'][index+1])            #add les 2 gènes au bloc en cours (set -> pas de doublons)

            else :                                              #si on est à la fin d'un cluster 
                if len(cluster) > 3 :                           #au moins 4 gènes     
                    tup=cluster,colName[0]                      #créer tuple set(cluster),E/C
                    E_strains.append(tup)                       #on l'ajoute à la liste de cluster asso à E/C selon colName
                cluster = set()                                 #et on reset le cluster en cours qu'il ait été add ou pas
                    
    if len(cluster) > 3 :                                       #si mon dernier cluster > 3 gènes
        tup=cluster,colName[0]                        
        E_strains.append(tup)                                   #on l'ajoute aussi (avec colName[0])

sorted_E_strains = sorted(E_strains,key=lambda x : len(x[0]),reverse=True)

#je veux [ [{cluster} , nbE , nbC] ] -> list of list :)

deduplicated_E_strains = []

for i in range(len(sorted_E_strains)):
    block,strain = sorted_E_strains[i]   
    index = [x for x,y in enumerate(deduplicated_E_strains) if y[0]==block]
    
    if len(index)==0 :
        if strain == 'E':
            deduplicated_E_strains.append([block,1,0])
        elif strain == 'C': 
            deduplicated_E_strains.append([block,0,1])
    
    elif len(index)==1:
        index = index[0]
        group = deduplicated_E_strains[index]
        
        if strain == 'E':
            group[1]=group[1]+1

        elif strain == 'C':
            group[2]=group[2]+1

    else :
        print('len index > 1 :', print(len(index)))
        
#sort par occ du 2eme element == E, puis par len si same occ
sorted_deduplicated_E_strains = sorted(deduplicated_E_strains,key=lambda x:(x[1],len(x[0])),reverse=True)

E_to_group = []

for block,E,C in sorted_deduplicated_E_strains :
    if E >= 2:
        E_to_group.append((block,E,C))


grouped_E_strains = []

for i in range(len(E_to_group)):
#    print('\n')
    
    block,E,C = E_to_group[i]
#    print('block',i+1,':',block)
    
    common = [element for element,E,C in grouped_E_strains if not block.isdisjoint(element)]
    
    if len(common)==0 :
        grouped_E_strains.append((block,E,C))
#        print('block',i+1,'( len',len(block),', occ E/C',E,C,') of E_to_group added to grouped list as block',grouped_E_strains.index((block,E,C))+1)

        
    elif len(common)==1:
        issubset = [element for element,E,C in grouped_E_strains if block.issubset(element)]
        issuperset = [element for element,E,C in grouped_E_strains if block.issuperset(element)]
        
        if len(issubset)==1:
            clust = issubset[0]
            clust_to_search = [tup for tup in grouped_E_strains if tup[0]==clust]
            block2,E2,C2 = clust_to_search[0]
#            print('block',i+1,'( len',len(block),', occ E/C',E,C,') of E_to_group is sub set of block', 
#                  grouped_E_strains.index((block2,E2,C2))+1,'( len',len(block2),', occ E/C',E2,C2,') in grouped list')
       
#        elif len(issubset)>1:
#            print('block is subset of',len(issubset),'other groups in grouped list:\n',issubset)
        
           
        if len(issuperset)==1:
            clust = issuperset[0]
            clust_to_search = [tup for tup in grouped_E_strains if tup[0]==clust]
            block2,E2,C2 = clust_to_search[0]
#            print('block',i+1,'( len',len(block),', occ E/C',E,C,') of C_to_group is super set of block', 
#                  grouped_E_strains.index((block2,E2,C2))+1,'( len',len(block2),', occ E/C',E2,C2,') in grouped list')
        
#        elif len(issuperset)>1:
#            print('block is superset of',len(issuperset),'other groups in grouped list:\n',issuperset)
                
    else : 
#        print('\ngenes in block',i+1,'( len',len(block),'occ E/C',E,C,') commons with',len(common),'other group.s :')
        for clust in common :
            clust_to_search = [tup for tup in grouped_E_strains if tup[0]==clust]
            block2,E2,C2 = clust_to_search[0]
#            print('\tblock',grouped_E_strains.index((block2,E2,C2))+1,'( len',len(block2),', occ E/C',E2,C2,') in grouped list, nb genes :',len(block.intersection(block2)))
        



    ########## generate a scatter plot to visualize the frequency of the enriched clusters ##########


gpa['Func_COG'] = gpa['Func_COG'].replace('-','S')

to_scatter={}
for cluster,occE,occC in grouped_C_strains :
    #get presence of each gene from cluster in E and C
    pres_C = []
    pres_E = []
    for gene in cluster : 
        line =  df_all.loc[df_all['Gene'] == gene]
        pres_E.append(line['pe'].values[0])
        pres_C.append(line['pc'].values[0])
    
    E = (sum(pres_E)/len(pres_E))/79*100
    C = (sum(pres_C)/len(pres_C))/64*100
    
    #cluster's distribution plot
    to_scatter['C'+str(grouped_C_strains.index((cluster,occE,occC))+1)]=(E,C,len(cluster))

    
for cluster,occE,occC in grouped_E_strains : 
    #get presence of each gene from cluster in E and C
    pres_C = []
    pres_E = []
    for gene in cluster : 
        line =  df_all.loc[df_all['Gene'] == gene]
        pres_E.append(line['pe'].values[0])
        pres_C.append(line['pc'].values[0])
    
    E = (sum(pres_E)/len(pres_E))/79*100
    C = (sum(pres_C)/len(pres_C))/64*100
    #cluster's distribution plot
    to_scatter['E'+str(grouped_E_strains.index((cluster,occE,occC))+1)]=(E,C,len(cluster))

fig, ax = plt.subplots(figsize=(9,7))


for k,(E,C,s) in to_scatter.items():
    if k.startswith('E'):
        c = 'green'
    else : 
        c = 'purple'
    sc = ax.scatter(E,C,s=s*5,label=k, c=c,alpha = 0.5)
    
plt.ylabel('% Presence in clinical isolates', fontsize=14)
plt.xlabel('% Presence in environmental isolates', fontsize=14)
plt.title('Frequencies of enriched clusters',fontsize=16)


green_patch = mpatches.Patch(color='green',alpha=0.5, label='Clusters from genes enriched in environmental isolates')
purple_patch = mpatches.Patch(color='purple',alpha=0.5, label='Clusters from genes enriched in clinical isolates')
handles=[green_patch,purple_patch]
ax.legend(handles=[green_patch,purple_patch])
