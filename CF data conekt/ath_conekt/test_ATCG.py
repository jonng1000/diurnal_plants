# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 16:00:23 2018

@author: jonathan

Tests to see if each DNA sequence only contains ATCG and if so, writes them
to a new fasta file 
"""

from Bio import Entrez, SeqIO

original_file = "ath_conekt.fa"
corrected_file = "ath_conekt_edited.fa"
ATCG = {'A', 'T', 'C', 'G'}
problems =[] 
with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        if set(record.seq) <= ATCG:
            SeqIO.write(record, corrected, 'fasta')
        else:
            problems.append(record)
            
# Second part, based upon results from 
problems_dict = {}
for one in problems:
    problems_dict[one.id] = set(one.seq)
