# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 23:56:47 2018

@author: 
    
This script extracts gene descriptions from Mercator and creates a file with gene IDs and
descriptions.
"""

import csv

original_file = "mercator.results.txt"
gene_desc_file = "name_description.txt"

with open(original_file) as original, open(gene_desc_file , 'w',  newline='') as file3:
    csvreader = csv.reader(original, delimiter='\t')
    csvwriter = csv.writer(file3, delimiter='\t')
    next(csvreader)
    for row in csvreader:
        if row[2].strip("'"):
            orig_gene = row[2].strip("'")
            gene_ID = orig_gene[0].upper() + orig_gene[1:]
            csvwriter.writerow([gene_ID, row[3].strip("'")])