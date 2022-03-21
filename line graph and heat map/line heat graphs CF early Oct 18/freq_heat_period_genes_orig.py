# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 15:44:05 2018

@author: weixiong
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import collections
import os
import csv
from scipy.stats import zscore

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
            # counts the frequencies of genes peaking at each hour,
            # sorts them according to the hour variable,
            count_h = collections.Counter(lower_than_p['LAG'])
            all_items = [(k,v) for k,v in count_h.items()]
            sorted_all = sorted(all_items, key=lambda one: one[0])
            hour = [x[0] for x in sorted_all]
            freq = [x[1] for x in sorted_all]
            # uses the hour and freq variables to plot a line graph
            plt.plot(hour, freq, '.-')
            plt.title(name)
            plt.savefig('line_' + name + '.pdf')
            plt.close()
            # z-scores table and plots a heatmap
            criteria = gene_ex_df['ADJ.P'] < 0.05
            lower_than_p = gene_ex_df[criteria]
            name = file.split('.')[1]
            sorted_lag = lower_than_p.sort_values(by=['LAG'])
            cleaned = sorted_lag.drop(['BH.Q', 'ADJ.P', 'PER', 'AMP', 'LAG'], axis=1)
            zcleaned = zscore(cleaned, axis=1)
            pic = sns.heatmap(zcleaned, cmap="bwr", xticklabels=False, yticklabels=False)
            plt.savefig('heatmap_' + name + '.pdf')
            plt.close()
