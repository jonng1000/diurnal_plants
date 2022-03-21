# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 00:27:48 2018

@author: weixiong

This script is to do some basic checking to see if there are problems with
the downloaded cds files and CF's JTK file.
"""

import pandas as pd
from Bio import SeqIO

# Checks to see if CF's JTK has only unique genes
df = pd.read_csv('JTK.Osa_Mat.txt', delimiter="\t", index_col=0)
all_idx = df.index.values # 42 330 genes

test_list = [x[:3] for x in all_idx]
test_set = set(test_list) # {'LOC', 'UR1', 'Uni'}
nonLOC_list = [x for x in all_idx if x[:3] != 'LOC'] # 3221 genes with 'UR1' or 'Uni'
# therefore, 42 330 - 3221 = #39 109 genes starting with 'LOC'
LOC_list = [x for x in all_idx if x[:3] == 'LOC'] # 39 109 genes
len(LOC_list) + len(UR1_list) + len(Uni_list) # 42 330 genes => same as all_idx

temp_list1 = [x.split('.')[0] for x in LOC_list]
temp_list2 = [x.split('.')[0] + '.' + x.split('.')[1] for x in UR1_list]
temp_set = set(temp_list1 + temp_list2 + Uni_list)
len(temp_set) #40 634 => CF JTK does not have all unique genes
len(set(LOC_list)) #39 109 => genes with 'LOC' are all unique

# Checks to see if CF's genes have the same IDs as downloaded cds files
record_dict1 = SeqIO.to_dict(SeqIO.parse("all.cds.fa", "fasta"))
keys_cds_LOC = [x for x in record_dict1.keys() if x[:3] != 'LOC']
keys_cds_nonLOC = [x.split('.')[0] for x in record_dict1.keys() if x[:3] == 'LOC']
cds_set = set(keys_cds_LOC + keys_cds_nonLOC)
len(cds_set) #55 986 genes

cds_inter_CF = cds_set & temp_set
len(cds_inter_CF) # 36 909 genes
len(set(temp_list1) - cds_set) # 513 genes which are in CF's LOC genes and are not in cds
