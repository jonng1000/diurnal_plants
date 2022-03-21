# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 23:56:47 2018

@author: weixiong

Edits kfl GO, InterPro and nane description file by adding _v1.1 back to the gene IDs.
"""

import csv
import pandas as pd

original_file = "GO_interproscan_pred.tsv"
edited_file = "GO_interproscan_pred_edit.tsv"

with open(original_file, newline='') as original, open(edited_file, 'w', newline='') as file3:
    csvreader = csv.reader(original, delimiter='\t')
    csvwriter = csv.writer(file3, delimiter='\t')
    for row in csvreader:
        gene_ID = row[0].lower() + '_v1.1'
        new_row = [gene_ID] + row[1:]
        csvwriter.writerow(new_row)
