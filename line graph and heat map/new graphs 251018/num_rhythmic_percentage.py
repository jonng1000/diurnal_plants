# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 15:44:05 2018

@author: weixiong
This script only finds out the number of rhythmic genes, total genes and
% rhythmic genes
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import collections
import os
import csv
from scipy.stats import zscore

with open('num_percent_rhythmic.txt', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['species', 'rhythmic genes', 'total genes', '% rhythmic'])
    # calculates % of rhythmic genes
    for file in os.listdir():
        if file.startswith('JTK'):
            gene_ex_df = pd.read_csv(file, delimiter="\t", index_col=0)
            criteria = gene_ex_df['ADJ.P'] < 0.05
            lower_than_p = gene_ex_df[criteria]
            sign_num = len(lower_than_p)/len(gene_ex_df)*100
            # writes sample name and % rhythmic genes to file
            name = file.split('.')[1]
            writer.writerow([name, len(lower_than_p), len(gene_ex_df), str(sign_num) + '%'])

