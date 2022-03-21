# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 16:38:57 2019

@author: jonathan
"""

###############################################################################
# From cds fasta file downloaded from the website, replaces headers with
# gene IDs (without decimal points), for uploading gene sequences to conekt 
# modified from original create_conekt_files.py
# us this on both peptide and cds fasta files
###############################################################################
from Bio import SeqIO

original_file = "all.cds.fa"

with open(original_file) as orig:
    records = SeqIO.parse(orig, 'fasta')
    csvwriter = csv.writer(new_desc, delimiter='\t')
    for record in records:
        break