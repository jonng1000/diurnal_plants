# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 11:47:45 2018

@author: weixiong

Generate gene clusters, which are defined by genes with the same lag values.
Sorts the cluster in ascending lag, where the first cluster is lag 0, second
has a higher lag and so on.

Outputs file called gene_clusters.txt
"""

input_f = 'JTK.DRA006158_edit.txt'
output_f = 'gene_clusters.txt'

import numpy as np
import pandas as pd
import csv

gene_ex_df = pd.read_csv(input_f, delimiter="\t", index_col=0)
criteria = gene_ex_df['ADJ.P'] < 0.05
lower_than_p = gene_ex_df[criteria]

LAG_set = set(gene_ex_df['LAG'])
LAG_dict = {}

for period in LAG_set:
    selection =lower_than_p['LAG'] == period
    selected =lower_than_p[selection]
    selected_genes = [item.split('.')[0] for item in selected.index]
    LAG_dict[period] = selected_genes
    
LAG_dkeys = list(LAG_dict.keys())
LAG_dkeys.sort()

with open(output_f, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    for k in range(24): #LAG_dkeys # also put a blank list below
        if k in LAG_dict:
            writer.writerow(LAG_dict[k])
        else:
            writer.writerow(['abc'])
