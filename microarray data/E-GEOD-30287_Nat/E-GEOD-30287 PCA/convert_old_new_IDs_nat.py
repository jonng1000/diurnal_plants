# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 12:40:10 2019

@author: weixiong
"""

import csv
import pandas as pd

compare_file = 'results_tab.txt'
gene_ex_file =  'E-GEOD-30287_combined_PCA.txt'
corrected_file = 'E-GEOD-30287_combined_PCA_edit.txt'

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
    # changes x to str(x)
    drop = [x for x in old_IDs if str(x) not in geneIDs_dict.keys()]
    corrected = gene_ex_df.drop(drop)
    
    cor_index = corrected.index.values
    # changes x to str(x)
    new_IDs = [geneIDs_dict[str(x)] for x in cor_index]
    corrected.index = new_IDs
    nr_genes = corrected[~corrected.index.duplicated(keep='first')]
    

nr_genes.to_csv(corrected_file, sep='\t')