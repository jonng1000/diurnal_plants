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

df = pd.read_csv('Smoellendorffii_91_v1.0.annotation_info.txt', delimiter="\t", index_col=0)
gene_desc_file = "name_description.txt"
records = SeqIO.parse('Smoellendorffii_91_v1.0.cds_primaryTranscriptOnly.fa', 'fasta')

record_dict = {record.description.split(' ')[1].split('=')[1]:record.id for record in records}

with open(gene_desc_file, 'w', newline='') as new_desc:
    csvwriter = csv.writer(new_desc, delimiter='\t')
    for pacID in df.index:
        name = record_dict[str(pacID)]
        if pd.isna(df.loc[pacID, 'arabi-defline']):
            desc = 'Function unknown'
        else:
            desc = df.loc[pacID, 'arabi-defline']
        csvwriter.writerow([name, desc])
