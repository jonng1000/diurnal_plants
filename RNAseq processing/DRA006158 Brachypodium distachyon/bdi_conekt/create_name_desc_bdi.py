# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 21:50:00 2019

@author: weixiong

Creates name_description.txt which is a file containing a list of gene IDs
and their description, from JGI Phytozome's annotation_info.txt file

Modified from create_name_desc
"""

import pandas as pd
from Bio import Entrez, SeqIO
import csv

df = pd.read_csv('Bdistachyon_314_v3.1.annotation_info.txt', delimiter="\t", index_col=0)
# Using original fasta file, not my edited one (removed decimals) as that will return
# many rows of descriptions corresponding to the same gene ID, and I only want one
# row of description, which needs me to use the transcript ID from the original
# file.
record_dict = SeqIO.to_dict(SeqIO.parse("Bdistachyon_314_v3.1.cds_primaryTranscriptOnly.fa", "fasta"))
gene_desc_file = "name_description.txt"

with open(gene_desc_file, 'w', newline='') as new_desc:
    csvwriter = csv.writer(new_desc, delimiter='\t')
    for key in record_dict.keys():
        temp = df.loc[df['transcriptName'] == key]
        # assumes a specific string format, need to change to suit specific files
        if key.count('.') == 1:
            name = key.split('.')[0]
        elif key.count('.') == 2:
            name = key.split('.')[0] + '.' + key.split('.')[1]
        else:
            raise ValueError('Unexpected number of dots')
        
        # checks if this is blank
        if not temp['arabi-defline'].any():
            desc = 'Function unknown'
        # selects for a string, as this returns an array
        else:
            desc = temp['arabi-defline'].values[0]
        csvwriter.writerow([name, desc])
