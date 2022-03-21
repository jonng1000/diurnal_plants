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

file = 'JTK.zea_may_edited.txt'
gene_ex_df = pd.read_csv(file, delimiter="\t", index_col=0)

criteria = gene_ex_df['ADJ.P'] < 0.05
lower_than_p = gene_ex_df[criteria]
name = file.split('.')[1]

sorted_lag = lower_than_p.sort_values(by=['LAG'])
cleaned = sorted_lag.drop(['BH.Q', 'ADJ.P', 'PER', 'AMP', 'LAG'], axis=1)
zcleaned = zscore(cleaned, axis=1)
zcleaned[zcleaned > 3] = 3
zcleaned[zcleaned < -3] = -3
pic = sns.heatmap(zcleaned, cmap="bwr", xticklabels=False, yticklabels=False)
#plt.savefig('heatmap_' + name + '.pdf')
plt.close()
sorted_lag.to_csv('zea_may_raw_score.txt', sep='\t')
np.savetxt('zea_may_zscore.txt', zcleaned, delimiter='\t')
