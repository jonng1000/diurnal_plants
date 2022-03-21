# -*- coding: utf-8 -*-
"""
Created on Fri May 25 18:13:07 2018

@author: weixiong001
"""
import os
import csv
import pandas as pd

for file in os.listdir('./'):
    if file.endswith('.txt'):
        input_file = file
        sample_name = input_file[:3]
        gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
        gene_ex_df_col_names = gene_ex_df.columns
        transposed = gene_ex_df_col_names.T
        
        with open('timepoints_in_cols.tsv', 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for line in transposed:
                writer.writerow([sample_name] + [line])