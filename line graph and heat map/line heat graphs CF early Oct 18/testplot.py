# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 11:17:32 2018

@author: weixiong
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import zscore

file = 'JTK.Syn_Mat.txt'
gene_ex_df = pd.read_csv(file, delimiter="\t", index_col=0)

criteria = gene_ex_df['ADJ.P'] < 0.05
lower_than_p = gene_ex_df[criteria]
name = file.split('.')[1]

sorted_lag = lower_than_p.sort_values(by=['LAG'])
cleaned = sorted_lag.drop(['BH.Q', 'ADJ.P', 'PER', 'AMP', 'LAG'], axis=1)
zcleaned = zscore(cleaned, axis=1)

pic = sns.heatmap(zcleaned, cmap="bwr", xticklabels=False, yticklabels=False)
#plt.savefig('heatmap_' + name + '.pdf')
plt.close()
#zcleaned.to_csv('Syn_Mat_zscore.txt', sep='\t')
#ran this: zcleaned.max()


'''
test = zscore(cleaned)

# supposed to work, apply zscore along rows
test2 = cleaned.apply(zscore, axis=1)
-> but get series object? error when generating heatmap

#converting to dataframe doesnt help
test2.to_frame()


#apply zscore along rows
test3 = zscore(cleaned, axis=1)
-> works, why is it different from apply function above?
gives numpy.ndarray datatype, in contrast to test2, which gives pandas.series datatype
'''