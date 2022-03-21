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
df = pd.read_csv('JTK.Kfl_Mat_edit.txt', delimiter="\t", index_col=0)
all_idx = df.index.values
temp_list = [x.split('_v1.1')[0] for x in all_idx]
temp_set = set(temp_list)
len(temp_set) #17053, same as len(temp_list) => CF JTK only has unique genes

record_dict1 = SeqIO.to_dict(SeqIO.parse("160614_klebsormidium_v1.1_CDS.fasta", "fasta"))
keys_cds = [x.split('_v1.1')[0] for x in record_dict1.keys()]
cds_set = set(keys_cds)
len(cds_set) #17283, same as len(keys_cds) => cds file only has unique genes

