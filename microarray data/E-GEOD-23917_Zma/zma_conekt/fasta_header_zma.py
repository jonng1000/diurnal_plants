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

#original_file = "Zmays_284_Ensembl-18_2010-01-MaizeSequence.protein_primaryTranscriptOnly.fa"
#corrected_file = "zma_aa_conekt.fa"
original_file = "Zmays_284_Ensembl-18_2010-01-MaizeSequence.cds_primaryTranscriptOnly.fa"
corrected_file = "zma_conekt.fa"

with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = record.description.split('locus=')[1].split(' ')[0]
        record.id = new_id
        record.description = new_id
        SeqIO.write(record, corrected, 'fasta')
        

