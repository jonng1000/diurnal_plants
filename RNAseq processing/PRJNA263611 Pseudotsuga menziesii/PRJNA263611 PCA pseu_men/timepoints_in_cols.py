# -*- coding: utf-8 -*-
"""
Created on Fri May 25 18:13:07 2018

@author: weixiong001
"""
'''
This is to select all time points in my gene expression matrix and transpose
them vertically
'''
import csv
import pandas as pd

input_file = 'PRJNA263611_combined.txt'

gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
gene_ex_df_col_names = gene_ex_df.columns
transposed = gene_ex_df_col_names.T

with open('timepoints_as_cols.tsv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    for line in transposed:
        writer.writerow([line])