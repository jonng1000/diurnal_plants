# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:54:29 2018

@author: weixiong

Creates a tab-delimited file with gene IDs (no decimals) in one column, and
gene description in another column
Uses the Biopython module
For uploading gene names and descriptions to conekt
No empty description for nat genes
"""
import csv
from Bio import SeqIO

original_file = "Nicotiana_attenuata.NIATTr2.cds.all.fa"
new_file = "name_description.txt"

with open(original_file) as original, open(new_file, 'w', newline='') as new:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        name = record.description.split('gene:')[1].split(' ')[0]
        desc = record.description.split('description:')[1]
        # Checks if there's any genes without a description
        if not desc:
            raise ValueError('Empty string')
        csvwriter = csv.writer(new, delimiter='\t')
        csvwriter.writerow([name, desc])
