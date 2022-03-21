# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 17:31:22 2018

@author: jonathan

This edits CF's ppu aa file so that gene IDs are in the same format as her
JTK output.
"""

from Bio import Entrez, SeqIO

original_file = "Porphyridium_genemodels_UPDATED.fasta"
corrected_file = "ppu_aa_conekt.fa"

with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        record.seq = record.seq.strip('*')
        SeqIO.write(record, corrected, 'fasta')
            