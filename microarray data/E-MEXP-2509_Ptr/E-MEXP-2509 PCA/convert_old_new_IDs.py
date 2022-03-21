# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 12:40:10 2019

@author: weixiong
"""

import csv
import pandas as pd

compare_file = 'results_tab.txt'
gene_ex_file =  'E-MEXP-2509_controls_PCA.txt'
corrected_file = 'E-MEXP-2509_controls_PCA_edit.txt'

# Creates key:value pairs as old gene ID: new gene ID pairs
with open(compare_file, newline='') as cfile:
    reader = csv.reader(cfile, delimiter='\t')
    geneIDs_dict = {}
    for line in reader:
        if line[0] not in geneIDs_dict:
            geneIDs_dict[line[0]] = line[1]
        else:
            continue

# Loads and replaces index in dataframe with correct gene IDs      
with open(gene_ex_file, newline='') as gefile:
    gene_ex_df = pd.read_csv(gefile, delimiter="\t", index_col=0)
    old_IDs = gene_ex_df.index.values
    drop = [x for x in old_IDs if x not in geneIDs_dict.keys()]
    corrected = gene_ex_df.drop(drop)
    
    cor_index = corrected.index.values
    new_IDs = [geneIDs_dict[x] for x in cor_index]
    corrected.index = new_IDs
    nr_genes = corrected[~corrected.index.duplicated(keep='first')]
    

nr_genes.to_csv(corrected_file, sep='\t')
# misc work to see why there are gene IDs in blast output (geneIDs_dict) which
# do not appear in the gene expression matrix (corrected). This means that
# the array probe file has sequences which do not appear in the gene expression
# matrix even though they map to actual zma genes. About 400 sequences fall into
# this category
    
# There are a total of ~ 30 000 probes which map to gene IDs in the gene
# expression matrix. Originally, the matrix has ~ 100 000 probes, which means
# that ~ 70 000 of them to not map to any gene IDs from BLAST.
    
# don't have to run this in future
'''
cor_set = set(corrected.index.values) # 30275
dict_set = set(geneIDs_dict.keys()) # 30680
diff = list(dict_set - cor_set) # 405
'''
# Checking how many genes are mapped by the probes, and how many of them are
# unique
'''
genes = [x for x in geneIDs_dict.values()] # 30680
len(set(genes)) # 21058

repeats = []
for x in genes:
    if genes.count(x) > 1:
        repeats.append(x)
'''