# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:00:54 2018

@author: weixiong
From cds fasta file downloaded from the website, replaces headers with gene IDs (without decimal points), for
uploading gene sequences to conekt 
"""

from Bio import SeqIO

original_file = "Solanum_lycopersicum.SL2.50.pep.all.fa"
corrected_file = "sly_aa_conekt.fa"

with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = record.id.split('.')[0]
        record.id = new_id
        record.description = new_id
        SeqIO.write(record, corrected, 'fasta')
        

