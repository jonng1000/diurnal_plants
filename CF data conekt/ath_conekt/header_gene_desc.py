# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 17:31:22 2018

@author: jonathan

This removes isoforms from genes, as it only selects transcripts with .1. Also
creates a file containing only gene name and its description.
"""

from Bio import Entrez, SeqIO
import csv

original_file = "Arabidopsis_thaliana.TAIR10.cds.all.fa"
corrected_file = "ath_conekt.fa"
gene_desc_file = "name_description.txt"

with open(original_file) as original, open(corrected_file, 'w') as corrected, open(gene_desc_file , 'w') as file3:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        if record.id.split('.')[1] == '1':
            try:
                gene_desc = record.description.split('description:')[1]
            except IndexError:
                gene_desc = 'Function unknown'
            
            temp_str = record.id.split('.')[0].lower()
            new_id = temp_str[0].upper() + temp_str[1:]
            record.id = new_id
            record.description = new_id
            SeqIO.write(record, corrected, 'fasta')
            
            csvwriter = csv.writer(file3, delimiter='\t')
            csvwriter.writerow([record.id, gene_desc])
            
            
            
            