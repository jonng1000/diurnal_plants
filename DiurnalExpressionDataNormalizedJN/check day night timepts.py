# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:07:23 2018

@author: weixiong001
"""
'''
This script is to check that all day/night timepoints make sense
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

for file in os.listdir('./'):
    if file.endswith('.txt'):
        input_file = file
        gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
        gene_ex_df_row_names = gene_ex_df.index
        gene_ex_df_col_names = gene_ex_df.columns

        list_timepoints = []
        for i in gene_ex_df_col_names:
            day_timepoint = i.split(" ")
            timepoint = day_timepoint[0]
            #processed_timepoint = timepoint.split(".")
            list_timepoints.append(timepoint)
        with open('all_timepoints.tsv', 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([input_file] + list_timepoints)
    