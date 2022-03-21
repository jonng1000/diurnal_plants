# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 17:31:22 2018

@author: jonathan

This ensures that gene ID is in lowercase, except for the first letter.
"""

from Bio import Entrez, SeqIO
import csv

original_file = "ath_aa_conekt.fa"
corrected_file = "ath_aa_conekt_edited.fa"

with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = record.id[0] + record.id[1:].lower()
        record.id = new_id
        record.description = new_id
        SeqIO.write(record, corrected, 'fasta')
