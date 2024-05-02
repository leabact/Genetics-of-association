#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 10:27:40 2022
@author: lmasson

comparison of genome length :
    generate all genome size (genome same = assembly file name)
    write size in list env or clin
    calculate GC%
    graph : boxplot
"""

from pathlib import Path as path
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as spicy

assemblies = path('/home/lmasson/Documents/stage_M2/genetics_of_association/assemblage_ref_ssg/assembly')

genome_size_E = []
genome_size_C = []
GC_E = []
GC_C = []


for f in path.iterdir(assemblies):
    genome_size = 0
    for seq_record in SeqIO.parse(f,'fasta'):
        genome_size = genome_size + len(seq_record.seq)
    if f.name.startswith('E'):
        genome_size_E.append(genome_size)
        GC_E.append(GC(seq_record.seq))
    else :
        genome_size_C.append(genome_size)
        GC_C.append(GC(seq_record.seq))


#size plot
size_to_plot = [genome_size_E,genome_size_C]

plt.figure(figsize=(5,4))

bp_size = plt.boxplot(size_to_plot,labels=['Environmental\nisolates', 'Clinical\nisolates'], patch_artist=True)
plt.ylabel('Genome size (million base pair)')
plt.title('Genome size',fontsize=14)
colors=['#CCFFE5','#CCCCFF']
for patch, color in zip(bp_size['boxes'], colors):
    patch.set_facecolor(color)
for median in bp_size['medians']:
    median.set(color ='black', linewidth = 1)

size_E = np.array(genome_size_E)
size_C = np.array(genome_size_C)

stat_size = spicy.ttest_ind(size_E,size_C)
pval = float(stat_size.pvalue)
if pval>0.05 :
    pval = 'ns'
elif 0.001 < pval <= 0.05 :
    pval='*'
elif pval <= 0.001 :
    pval = '***'

plt.tick_params(axis='both', labelsize = 12)

plt.hlines(y=7380000,xmin =1, xmax=2, linewidth=1, color ='grey')
plt.text(1.5 ,7300000 ,pval , color = 'black', fontvariant = 'small-caps', ha = 'center', fontsize = 12) 
plt.show(bp_size)

'''
#GC plot + stat
GC_to_plot = [GC_E,GC_C]
bp_gc = plt.boxplot(GC_to_plot,labels=['Environmental\nsamples', 'Clinical\nsamples'], patch_artist=True)
plt.ylabel('GC%')
plt.title('GC%')
colors=['#CCFFE5','#CCCCFF']
for patch, color in zip(bp_gc['boxes'], colors):
    patch.set_facecolor(color)
for median in bp_gc['medians']:
    median.set(color ='black', linewidth = 1)

np_GC_E = np.array(GC_E)
np_GC_C = np.array(GC_C)

stat_GC = spicy.ttest_ind(np_GC_E,np_GC_C)
pval = float(stat_GC.pvalue)

if pval>0.05 :
    pval = 'ns'
elif 0.001 < pval <= 0.05 :
    pval='*'
elif pval <= 0.001 :
    pval = '***'
    
plt.text(1.5, 75, pval, color = 'black', fontvariant = 'small-caps', ha = 'center', fontsize = 9)
plt.show(bp_gc)
'''

"""
res w/ stats: 
    mean_size_E : 6942194.94 -> 6.94 10⁶
    mean_size_C : 6601198.34 -> 6.60 10⁶
    mean_GC_E : 63.86
    mean_GC_C : 65.35
"""

