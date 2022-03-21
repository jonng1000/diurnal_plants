# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 21:50:00 2019

@author: weixiong

Creates name_description.txt which is a file containing a list of gene IDs
and their description, from JGI Phytozome's annotation_info.txt file.
This is named beta.py because a cleaner version is found in the
parent folder, named create_name_desc.py.
"""

import pandas as pd
from Bio import Entrez, SeqIO
import csv

# Checks to see if CF's JTK has only unique genes
df = pd.read_csv('Creinhardtii_281_v5.6.annotation_info.txt', delimiter="\t", index_col=0)
record_dict = SeqIO.to_dict(SeqIO.parse("Creinhardtii_281_v5.6.cds_primaryTranscriptOnly.fa", "fasta"))
gene_desc_file = "name_description.txt"

# Checks if these two columns are exactly the same
check = df['transcriptName'] ==	df['peptideName']
check.all() # True => these two columns are exactly the same

with open(gene_desc_file, 'w', newline='') as new_desc:
    csvwriter = csv.writer(new_desc, delimiter='\t')
    for key in record_dict.keys():
        temp = df.loc[df['transcriptName'] == key]
        name = key.split('.')[0] + '.' + key.split('.')[1]
        if not temp['arabi-defline'].any():
            desc = 'Function unknown'
        else:
            desc = temp['arabi-defline'].values[0]
        csvwriter.writerow([name, desc])
