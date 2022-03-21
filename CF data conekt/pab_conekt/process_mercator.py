# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 23:56:47 2018

@author: 
    
This script extracts gene descriptions from Mercator and creates a file with gene IDs and
descriptions. Modified from process_mercator.py to suit pab data.
"""

import csv

original_file = "pab_mercator.results.txt"
gene_desc_file = "name_description.txt"

with open(original_file) as original, open(gene_desc_file , 'w') as file3:
    csvreader = csv.reader(original, delimiter='\t')
    csvwriter = csv.writer(file3, delimiter='\t')
    next(csvreader)
    for row in csvreader:
        if row[2].strip("'"):
            gene_ID = row[2].strip("'")
            new = gene_ID[:2].upper() + gene_ID[2:]
            csvwriter.writerow([new, row[3].strip("'")])