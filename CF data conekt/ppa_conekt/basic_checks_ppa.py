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
df = pd.read_csv('JTK.Ppa_Mat_edit.txt', delimiter="\t", index_col=0)
all_idx = df.index.values # 32 275 genes

temp_list1 = [x.split('.')[0] for x in all_idx]
temp_set = set(temp_list1)
len(temp_set) #32 275 => CF JTK has all unique genes
