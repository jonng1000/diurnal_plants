# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 17:31:22 2018

@author: jonathan

This edits CF's Ppu cds file so that gene IDs are in the same format as her
JTK output.
"""

from Bio import Entrez, SeqIO

original_file = "Smoellendorffii_91_v1.0.protein_primaryTranscriptOnly.fa"
corrected_file = "smo_aa_conekt.fa"

with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = record.description.split(' ')[0]
        record.id = new_id
        record.description = new_id
        SeqIO.write(record, corrected, 'fasta')
            