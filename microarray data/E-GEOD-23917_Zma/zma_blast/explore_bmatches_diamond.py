# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 17:51:48 2019

@author: weixiong

Use to explore data from diamond.
See how many unique probes there are, and
how many are tied for lowest e-values.
"""

import pandas as pd
import csv
file = 'zma_matches.txt'

# 45 456 probes
df = pd.read_csv('zma_matches.txt', delimiter="\t", header=None, index_col=0)
index_list = df.index.tolist()
index_set = set(index_list)
len(index_set) # 22 224 unique probes

# 5 largest e-values are all 0.0009, values typically range around
# e-6
df.nlargest(5, 10)

probe_dict = {}
smaller_dict = []

with open(file, newline='') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        probe = row[0]
        e_v = float(row[10])
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
# ~41 000 probes and ~6 000 of them match to >1 gene with e-values tied for
# lowest value
# Diamond output is not as good as blastn, probably because its a more
# indirect way of comparing the nt probes to the protein database, as blastn
# directly compares the nt probes to the nt database