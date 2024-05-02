#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 09:40:17 2022
@author: lmasson

homemade python :
    parse gene_presence_absence : get the core genome + list of all genes
        establish core/softcore/shell/cloud gene groups : core gene (>=99%) (4117), 
        softcore gene (95-99%) (819), shell gene (15-95%) (2068), cloud gene (0-15%) (17512)
    calculate counting tabs (x = pres/abs, y = env/clin) -> don't take gene in env cloud+clin cloud / same with core
    calculate the pvalue with a ficher exact test + corrects it
    write results in all_res dataframe + separates gene enriched in the 2 pop° (df enriched_E / C)
    save res df with pickle : all_res, enriched_all, enriched_E, enriched_C
"""

import pandas as pd
import scipy.stats as spicy
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

############################# functions ##############################

def tab_modif (df) :
#    df.columns = df.columns.str.strip().str.replace(' ','_').str.replace('.','')
    df = df.fillna('no')
    return(df)


def core_genome() :
    core_genes = {}
    for line in range(len(df)):
        if df['No_isolates'][line] >= 142 : 
            core_genes[df['Gene'][line]]=df['Annotation'][line]
    return(core_genes)


def get_genes() :
    '''
    returns a list of gene in the defined group
    mini-maxi for gene group : softcore : 136-141, shell : 21-135, cloud : 0-20
    '''
    all_genes = []
    for index in range(len(df)) : 
            all_genes.append(df['Gene'][index])            
    return(all_genes)


def alexis_uwu(pe,pc,ae,ac):
    '''takes data set of presence/absence, returns it into df '''
    data = []
    data.append([pe, ae])
    data.append([pc,ac])
    df_count = pd.DataFrame(data, index=['Env','Clin'], columns=['pres','abs'])
    return(df_count)


def counting_tabs(all_genes):
    '''
    create sub df E/C from gene_presence_absence
    count presence/absence in env or clin
    from presence count : determine group of genes, separating env/clin
    return a counting table per gene, stocked in a dict {gene : df, gene : df, ...} 
    + lists of variables to write in df res
    '''
    select = [str(a).startswith("E") for a in df.columns]
    df_E = df[df.columns[select]]
    select = [str(a).startswith("C") for a in df.columns]
    df_C = df[df.columns[select]]

    tab_counts = {}
    gene_list = []
    annotation_list = []
    pe_list = []
    pc_list = []
    ae_list = []
    ac_list = []

    for gene in all_genes :
        select = df['Gene'] == gene
        if sum(select) != 1:
            print("warnings, pas le bon nbre de lignes: "+ gene)
        e = df_E.loc[select].squeeze()
        pe = sum([1 for i in e if i != 'no'])
        ae = len(df_E.columns) - pe
        c = df_C.loc[select].squeeze()
        pc = sum([1 for i in c if i != 'no'])
        ac = len(df_C.columns) - pc
        
        #if gene = cloud/core from E & C, skip it        
        if pe < 15 and pc < 10:
            continue
        if pe > 78 and pc > 63 :
            continue

        gene_list.append(gene)
        annotation_list.append(df.loc[df.Gene == gene, 'Annotation'].values[0])
        pe_list.append(pe)
        pc_list.append(pc)
        ae_list.append(ae)
        ac_list.append(ac)
        df_count = alexis_uwu(pe,pc,ae,ac)
        tab_counts[gene]=df_count

    return(tab_counts,gene_list,annotation_list,pe_list,pc_list,ae_list,ac_list)


def gene_enrichment_stats(tab_counts):
    '''
    do a fisher exact test on counting tabs
    correct the p values (Benjamini Hochberg procedure)
    return lists of pvalues and q values
    '''
    #calculates p values
    pval = []
    qval = []
    for gene,tab in tab_counts.items():
        a,b = spicy.fisher_exact(tab)
        pval.append(b)
    
    #calculate q values
    rejected, qvalues = fdrcorrection(pval)
    for qvalue in qvalues :
        qval.append(qvalue)
  
    return(pval,qval)


def to_df_res(gene_list, annotation_list, pe_list, pc_list, ae_list, ac_list, pval, qval):
    '''
    create df_res and add the data, returns full df_res
    '''
    df_res = pd.DataFrame(columns=['Gene','Annotation','pe','ae','pc','ac','pval','qval','enriched'])
    df_res['Gene'] = gene_list
    df_res['Annotation'] = annotation_list
    df_res['pe'] = pe_list
    df_res['ae'] = ae_list
    df_res['pc'] = pc_list
    df_res['ac'] = ac_list
    df_res['pval'] = pval
    df_res['qval'] = qval
    return(df_res)


def fill_enriched(df_res):
    '''
    if qval <= 0.05
    fill the enriched column in df res with :
        clin if pc > pe
        env if pe > pc
    '''
    for index in range(len(df_res)):
        if df_res['qval'][index] <= 0.05 :
            if df_res['pe'][index] > df_res['pc'][index] :
                df_res.loc[index,'enriched']='Env'                
            elif df_res['pc'][index] > df_res['pe'][index] :
                df_res.loc[index,'enriched']='Clin'
    df_res = df_res.fillna('no')
    
    return(df_res)


def get_enriched(df_res) : 
    '''
    create and returns 3 dataframes containing :
       df enriched : all genes statistically enriched in a environment
       df_E : genes enriched in environmental isolates
       df_C : genes enriched in clinical isolates 
    '''
    df_enriched = pd.DataFrame(columns=['Gene','Annotation','pe','pc','ae','ac','pval','qval','enriched'])
    df_E = pd.DataFrame(columns=['Gene','Annotation','pe','pc','ae','ac','pval','qval'])
    df_C = pd.DataFrame(columns=['Gene','Annotation','pe','pc','ae','ac','pval','qval'])
    
    for index in range(len(df_res)):
        if df_res['enriched'][index] != 'no':
            data = df_res.iloc[index]
            df_enriched = df_enriched.append(data, ignore_index = True)
            
            if df_res['enriched'][index] == 'Env' :
                 data = df_res.iloc[index]
                 df_E = df_E.append(data, ignore_index = True)
            if df_res['enriched'][index] == 'Clin' :
                 data = df_res.iloc[index]
                 df_C = df_C.append(data, ignore_index = True)

    return(df_enriched, df_E, df_C)

#add "group" col to df_E/df_C 
def write_gene_group(df, mini, maxi, core, pe_or_pc):
    '''
    create a new column in df, write the gene group in it
    env : mini = 15, maxi = 74, core = 78
    clin : mini = 10, maxi = 60, core = 63
    returns filled df
    '''
    for index in range(len(df)):
        if df[pe_or_pc][index] < mini :
            df.loc[index,'group']='cloud'    
            
        elif mini <= df[pe_or_pc][index] < maxi : 
            df.loc[index,'group']='shell' 
            
        elif maxi <= df[pe_or_pc][index] < core :
            df.loc[index,'group']='soft' 
  
        elif core <= df[pe_or_pc][index] :
            df.loc[index,'group']='core'

    return(df)


def gene_enriched_dict(df_enriched):
    '''
    separates genes enriched in clin or env
    returns dict {gene : annotation...}
    '''
    gene_enriched_env = {}
    gene_enriched_clin = {}
    for index in range(len(df_enriched)):
        if df_enriched['enriched'][index] == 'Env':
            gene_enriched_env[df_enriched['Gene'][index]] = df_enriched['Annotation'][index]
        elif df_enriched['enriched'][index] == 'Clin':
            gene_enriched_clin[df_enriched['Gene'][index]] = df_enriched['Annotation'][index]

    return(gene_enriched_env,gene_enriched_clin)


###plots###

def plotting_pqval(df_res, name):
    plt.hist(df_res['pval'],alpha=1, color = 'lightblue', bins = 20)
    plt.hist(df_res['qval'],alpha=0.4, color = 'orange', bins = 20)
    plt.axvline(x=0.05, color = 'grey', alpha = 0.5)
    plt.legend(['0.05','p values','q values'])
    plt.title(name+' : raw and corrected pvalues \nnumber of genes with qval <= 0.05 : '+str(sum(df_res['qval']<=0.05)))
    return(plt.show())


def scatter_all_genes(all_res):
    '''
    scatter plot all genes, with different color in enriched in clin, env or neither
    returns the figure
    '''
    res_number = all_res
    for index in range(len(res_number)):
        if res_number['enriched'][index] == 'no':
            res_number.loc[index,'enriched'] = 0
        elif res_number['enriched'][index] == 'Env':
            res_number.loc[index,'enriched'] = 1
        elif res_number['enriched'][index] == 'Clin':
            res_number.loc[index,'enriched'] = 2
    
    colors={0:'#B6FFFF',1:'#B9FFB9', 2:'#CCE5FF'}
    res_number.plot.scatter(x='pe',y='pc', c=res_number['enriched'].map(colors))
    no = mpatches.Patch(color='#B6FFFF', label = 'Non-enriched genes')
    clin = mpatches.Patch(color='#CCE5FF', label = 'Enriched in clinical')
    env = mpatches.Patch(color='#B9FFB9', label = 'Enriched in sink traps')
    plt.legend(handles=[no,env,clin], prop={'size':9},loc='upper left')
    plt.title('Genes repartition in isolates')
    plt.ylabel('Presence in clinical isolates')
    plt.xlabel('Presence in environmental isolates')
    
    return(plt.show())


def dflist_to_pickle(df_list):
    '''
    sereliaze df res to pickle
    load ex : enriched_all = pd.read_pickle('enriched_all.pickle')
    '''
    for df in df_list:
        name = [x for x in globals() if globals()[x] is df][0] 
        with open('pickle_stuff/'+name+'.pickle','wb') as f :
            df.to_pickle(f)
    return()


################################# MAIN ##################################

        #open+modify csv

#df = pd.read_csv('/home/lmasson/Documents/stage_M2/genetics_of_association/homemade_python/gene_presence_absence.csv', header = 0)
df = pd.read_csv('gene_infos/gpa.csv',header=0)
df = tab_modif(df)
print('modif tab ok')

        #get core genes and all genes
core_genes = core_genome()
all_genes = get_genes()

        #generate counting tabs and gene groups
#tabs :       pres   abs
#        env   pe    ae
#       clin   pc    ac
tab_counts, gene_list, annotation_list, pe_list, pc_list, ae_list, ac_list = counting_tabs(all_genes)
print('counting tabs ok')

        #stats for gene enrichment
pval, qval = gene_enrichment_stats(tab_counts)
print('enrichment calc ok')


        #write res in df, fill 'enriched' col with env/clin if qval <= 0.05
        #get enriched genes in new df
all_res = to_df_res(gene_list, annotation_list, pe_list, pc_list, ae_list, ac_list, pval, qval)
all_res = fill_enriched(all_res)
enriched_all, enriched_E, enriched_C = get_enriched(all_res)
print('df res ok')


        #in df enriched E/C, write the group of each gene
enriched_E = write_gene_group(enriched_E, 15, 74, 78, 'pe')
enriched_C = write_gene_group(enriched_C, 10, 60, 63, 'pc')


        #plot stuff
fig_scatter_all_genes = scatter_all_genes(all_res)
fig_pqval = plotting_pqval(all_res, 'all genes')

wedgeprops={"edgecolor":"black",'linewidth': 0.3, 'antialiased': True}

labelsenv = ['shell \n(880)','core (45)','soft (7)']
fig_group_repartition_E = enriched_E['group'].value_counts().plot(kind='pie', colormap = 'Pastel2', explode = (0,0.15,0),wedgeprops=wedgeprops,
                                              counterclock=False, textprops={'fontsize': 14}, labels=labelsenv)
fig_group_repartition_E.set_title('Environmentally enriched gene groups\n'+'Total enriched genes: '+str(len(enriched_E)))
fig_group_repartition_E.axis('equal')
plt.show(fig_group_repartition_E)

labelsclin = ['shell \n(130)','soft (2)']
fig_group_repartition_C = enriched_C['group'].value_counts().plot(kind='pie', colormap = 'Pastel2', wedgeprops=wedgeprops, counterclock=False,
                                              textprops={'fontsize': 14}, labels=labelsclin)
fig_group_repartition_C.set_title('Clinically enriched gene groups\n'+'Total enriched genes: '+str(len(enriched_C)))
fig_group_repartition_C.axis('equal')
plt.show(fig_group_repartition_C)

        #save df results with pickle
all_df_res = [all_res, enriched_all, enriched_E, enriched_C]
dflist_to_pickle(all_df_res)


