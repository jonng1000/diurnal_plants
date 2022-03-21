# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 23:56:47 2018

@author: weixiong

This script extracts gene descriptions from Mercator and creates a file with gene IDs and
descriptions.
"""

import csv
import pandas as pd

original_file = "ppu_mercator.results.txt"
gene_desc_file = "name_description.txt"

desc_dict = {}
with open(original_file, newline='') as original, open(gene_desc_file, 'w', newline='') as file3:
    csvreader = csv.reader(original, delimiter='\t')
    csvwriter = csv.writer(file3, delimiter='\t')
    next(csvreader)
    for row in csvreader:
        if row[2].strip("'") and row[2].strip("'").startswith('popu|evm'):
            gene_ID = row[2].strip("'")
            new = gene_ID.split('|')[1]
            desc = row[3].strip("'")
            if new not in desc_dict:
                desc_dict[new] = desc
                csvwriter.writerow([new, desc])
