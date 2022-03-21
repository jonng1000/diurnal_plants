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
df = pd.read_csv('JTK.Ppu_Mat_edit.txt', delimiter="\t", index_col=0)
all_idx = df.index.values # 8355 genes

temp_list1 = [x.split('.') for x in all_idx]
temp_list2 = ['.'.join(x[:3]) for x in temp_list1]
temp_set = set(temp_list2)
len(temp_set) #1531 => CF JTK does not have all unique genes, if decimals are removed
