# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 22:01:09 2019

@author: weixiong

Finds lengths of sequences in a fasta file
"""

import pandas as pd
from Bio import SeqIO

original_file = 'spe_conekt.fa'

with open(original_file) as orig:
    records = SeqIO.parse(original_file, 'fasta')
    length = 0
    len_list = []
    for record in records:
        if len(record.seq) > length:
            gene = record.id
            length = len(record.seq)
        len_list.append(len(record.seq))
        
len_list.sort(reverse=True)