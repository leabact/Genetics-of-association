#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
boxplot comparaison nombre de gènes env et clin
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as spicy

df = pd.read_csv('../gene_infos/gene_presence_absence.csv', header = 0)

df.columns = df.columns.str.strip().str.replace(' ','_').str.replace('.','')
df = df.drop(['Non-unique_Gene_name','No_sequences','Avg_sequences_per_isolate','Genome_Fragment',
              'Order_within_Fragment','Accessory_Fragment','Accessory_Order_with_Fragment','QC',
              'Min_group_size_nuc','Max_group_size_nuc','Avg_group_size_nuc'], axis = 1)

lst_genes_E = []
lst_genes_C = []

for (index,colName) in enumerate(df):
    if colName.startswith('E'):
        lst_genes_E.append(df[colName].notnull().sum())
    if colName.startswith('C'):
        lst_genes_C.append(df[colName].notnull().sum())

#nb gene plot
nb_genes_to_plot = [lst_genes_E, lst_genes_C]
plt.figure(figsize=(5,4))

bp_genes = plt.boxplot(nb_genes_to_plot,labels=['Environmental\nisolates', 'Clinical\nisolates'], patch_artist=True)
plt.ylabel('Number of genes')
plt.title('Number of genes per strain',fontsize=14)
colors=['#CCFFE5','#CCCCFF','#FFFFCC']
for patch, color in zip(bp_genes['boxes'], colors):
    patch.set_facecolor(color)
for median in bp_genes['medians']:
    median.set(color ='black', linewidth = 1)

#tab pour stat avec numpy
genes_E = np.array(lst_genes_E)
genes_C = np.array(lst_genes_C)

#stat avec spicy
stat_size = spicy.ttest_ind(genes_E,genes_C)
pval = float(stat_size.pvalue)
if pval>0.05 :
    pval = 'ns'
elif 0.001 < pval <= 0.05 :
    pval='*'
elif pval <= 0.001 :
    pval = '***'

plt.tick_params(axis='both', labelsize = 12)

plt.hlines(y=6900,xmin =1, xmax=2, linewidth=1, color ='grey')
plt.text(1.5 ,6800 ,pval, color = 'black', fontvariant = 'small-caps', ha = 'center', fontsize = 12) 
plt.show(bp_genes)

#what does it take to be a PA -> 5681 genes PAO1 (total gene cluster)
#ttest_ind(lst_genes_E, lst_genes_C) -> same res as spicy
