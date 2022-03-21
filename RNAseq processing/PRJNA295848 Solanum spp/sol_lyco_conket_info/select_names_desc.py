# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:54:29 2018

@author: weixiong

Creates a tab-delimited file with gene IDs (no decimals) in one column, and
gene description in another column
Uses the Biopython module
For uploading gene names and descriptions to conekt
"""
import csv
from Bio import SeqIO

original_file = "ITAG3.2_CDS.fasta"
new_file = "name_description.txt"

with open(original_file) as original, open(new_file, 'w', newline='') as new:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        name = record.id.split('.')[0]
        try:
            desc = record.description.split(' ', 1)[1]
        except IndexError:
            desc = 'Function unknown'

        csvwriter = csv.writer(new, delimiter='\t')
        csvwriter.writerow([name, desc])
