# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 15:44:05 2018

@author: weixiong

creates countplots, which is seaborn's name for bar graphs without
estimates of central tendency 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import csv

with open('percent_rhythmic.txt', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    # calculates % of rhythmic genes
    for file in os.listdir():
        if file.startswith('JTK'):
            gene_ex_df = pd.read_csv(file, delimiter="\t", index_col=0)
            criteria = gene_ex_df['ADJ.P'] < 0.05
            lower_than_p = gene_ex_df[criteria]
            sign_num = len(lower_than_p)/len(gene_ex_df)*100
            # writes sample name and % rhythmic genes to file
            name = file.split('.')[1]
            writer.writerow([name, str(sign_num) + '%'])
            
            sns.countplot(x='LAG', data=lower_than_p)
            plt.savefig('countplot_' + name + '.png')
            plt.close()
