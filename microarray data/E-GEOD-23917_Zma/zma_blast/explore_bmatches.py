# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 17:51:48 2019

@author: weixiong

Use to explore data from blast.
See how many unique probes there are, and
how many are tied for lowest e-values.
"""

import pandas as pd
import csv
file = 'results_tab.txt'

# 45 456 probes
df = pd.read_csv('results_tab.txt', delimiter="\t", header=None, index_col=0)
index_list = df.index.tolist()
index_set = set(index_list)
len(index_set) #30 680 unique probes

# 5 largest e-values are all 2.440000e-07, values typically range from
# e-10 to e-20
df.nlargest(5, 4)

probe_dict = {}
smaller_dict = []

with open(file, newline='') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        probe = row[0]
        e_v = float(row[4])
        if probe not in probe_dict:
            probe_dict[probe] = e_v
        else:
            if e_v > probe_dict[probe]:
                continue
            elif e_v == probe_dict[probe]:
                smaller_dict.append([probe, e_v])
            elif e_v < probe_dict[probe]:
                raise ValueError('Does not make sense')
                
# Conclusion
# ~45 000 probes and ~4 000 of them match to >1 gene with e-values tied for
# lowest value
# Discard these ~ 4000 matches