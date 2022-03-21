# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:00:54 2018

@author: weixiong
"""

###############################################################################
# From cds fasta file downloaded from the website, replaces headers with
# gene IDs (without decimal points), for uploading gene sequences to conekt 
# modified from original create_conekt_files.py
# us this on both peptide and cds fasta files
###############################################################################
from Bio import SeqIO

original_file = "Ptrichocarpa_444_v3.1.protein_primaryTranscriptOnly.fa"
corrected_file = "ptr_aa_conekt.fa"
#original_file = "Ptrichocarpa_444_v3.1.cds_primaryTranscriptOnly.fa"
#corrected_file = "ptr_conekt.fa"

with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = record.id.split('.')[0] + '.' + record.id.split('.')[1]
        record.id = new_id
        record.description = new_id
        SeqIO.write(record, corrected, 'fasta')
        
# Just to check if primary transcript file has only unique genes
# Yes it has
record_dict = SeqIO.to_dict(SeqIO.parse(corrected_file, "fasta"))
test = [key for key in record_dict.keys()] # 42950 gens
test_set = set(test)
len(test_set) # 42950 genes
