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
df = pd.read_csv('JTK.Cre_Mat.txt', delimiter="\t", index_col=0)
all_idx = df.index.values
temp_list = [x.split('.')[0] + '.' + x.split('.')[1] for x in all_idx]
temp_set = set(temp_list)
len(temp_set) #17733, same as len(temp_list) => CF JTK only has unique genes

# Checks to see if CF's genes have the same IDs as downloaded cds files ve r5.6
record_dict1 = SeqIO.to_dict(SeqIO.parse("Creinhardtii_281_v5.6.cds.fa", "fasta"))
record_dict2 = SeqIO.to_dict(SeqIO.parse("Creinhardtii_281_v5.6.cds_primaryTranscriptOnly.fa", "fasta"))

keys_cds = [x.split('.')[0] + '.' + x.split('.')[1] for x in record_dict1.keys()]
keys2_pri = [x.split('.')[0] + '.' + x.split('.')[1] for x in record_dict2.keys()]

cds_set = set(keys_cds)
len(cds_set) #17741, same as len(temp_list) => CF JTK only has unique genes

pri_set = set(keys2_pri)
len(pri_set) #17741, same as len(temp_list) => CF JTK only has unique genes

cds_inter_CF = cds_set & temp_set
len(cds_inter_CF) # 17733 genes
temp_set <= cds_inter_CF #True, CF genes is a subset of cds
pri_inter_CF = pri_set & temp_set
len(pri_inter_CF) # 17733 genes
temp_set <= pri_inter_CF #True, CF genes is a subset of pri

# Checks to see if CF's genes have the same IDs as downloaded cds files ve r5.5
record_dict1_55 = SeqIO.to_dict(SeqIO.parse("Creinhardtii_281_v5.5.cds.fa", "fasta"))
record_dict2_55 = SeqIO.to_dict(SeqIO.parse("Creinhardtii_281_v5.5.cds_primaryTranscriptOnly.fa", "fasta"))

keys_cds_55 = [x.split('.')[0] + '.' + x.split('.')[1] for x in record_dict1_55.keys()]
keys2_pri_55 = [x.split('.')[0] + '.' + x.split('.')[1] for x in record_dict2_55.keys()]

cds_set_55 = set(keys_cds_55)
len(cds_set_55) #17744, same as len(temp_list) => CF JTK only has unique genes

pri_set_55 = set(keys2_pri_55)
len(pri_set_55) #17741, same as len(temp_list) => CF JTK only has unique genes

cds55_inter_CF = cds_set_55 & temp_set
len(cds55_inter_CF) # 17733 genes
temp_set <= cds55_inter_CF #True, CF genes is a subset of cds
pri55_inter_CF = pri_set_55 & temp_set
len(pri_inter_CF) # 17733 genes
temp_set <= pri_inter_CF #True, CF genes is a subset of pri


diff_seq = [x for x in record_dict1.keys() if x not in record_dict2.keys()]
