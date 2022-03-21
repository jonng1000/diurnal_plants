# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 23:56:47 2018

@author: weixiong

This script extracts gene descriptions from Mercator and creates a file with gene IDs and
descriptions.
"""

import csv
import pandas as pd

original_file = "mercator.results.txt"
gene_desc_file = "name_description.txt"

desc_dict = {}
with open(original_file, newline='') as original, open(gene_desc_file, 'w', newline='') as file3:
    csvreader = csv.reader(original, delimiter='\t')
    csvwriter = csv.writer(file3, delimiter='\t')
    next(csvreader)
    for row in csvreader:
        if row[2].strip("'"):
            gene_ID = row[2].strip("'").upper()
            desc = row[3].strip("'")
            if gene_ID not in desc_dict:
                desc_dict[gene_ID] = desc
                csvwriter.writerow([gene_ID, desc])
