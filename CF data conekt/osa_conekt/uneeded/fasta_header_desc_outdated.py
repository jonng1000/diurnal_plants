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

# This script also generates the .txt file with gene names and descriptions
# Used for rice, which has inconsistent gene names in the cds file, hence
# script accounts for that

# Realised that its not good to combine both modifying fasta header and
# creating gene description in one script, as I use this for gene cds files,
# but don't need to create gene descriptions for protein fasta files.
###############################################################################
from Bio import SeqIO
import csv

original_file = "all.cds.fa"
corrected_file = "osa_conekt.fa"

gene_desc_file = 'name_description.txt'

with open(original_file) as orig, open(corrected_file, 'w') as corrected, open(gene_desc_file, 'w', newline='') as new_desc:
    records = SeqIO.parse(orig, 'fasta')
    csvwriter = csv.writer(new_desc, delimiter='\t')
    for record in records:
        # processes genes starting with LOC
        if record.id.startswith('LOC'):
            if record.id.split('.')[1] == '1':
                new_id = record.id.split('.')[0]
                record.id = new_id
                try:
                    description = record.description.split('|')[1]
                except IndexError:
                    description = 'Function unknown'
                record.description = new_id
                SeqIO.write(record, corrected, 'fasta')
                csvwriter.writerow([new_id, description])
        # processes genes starting with Chr
        elif record.id.startswith('Chr'):
            try:
                description = record.description.split('|')[1]
            except IndexError:
                description = 'Function unknown'
            record.description = record.id
            SeqIO.write(record, corrected, 'fasta')
            csvwriter.writerow([record.id, description])
        else:
            raise ValueError('Unknown gene name encountered')
            
###############################################################################
# Run above code block separately from this
# This code block modifies fasta headers for peptide files
###############################################################################

original_file = "all.pep"
corrected_file = "osa_aa_conekt.fa"
with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(orig, 'fasta')
    for record in records:
        # processes genes starting with LOC
        if record.id.startswith('LOC'):
            if record.id.split('.')[1] == '1':
                new_id = record.id.split('.')[0]
                record.id = new_id
                record.description = new_id
                SeqIO.write(record, corrected, 'fasta')
        # processes genes starting with Chr
        elif record.id.startswith('Chr'):
            record.description = record.id
            SeqIO.write(record, corrected, 'fasta')
        else:
            raise ValueError('Unknown gene name encountered')