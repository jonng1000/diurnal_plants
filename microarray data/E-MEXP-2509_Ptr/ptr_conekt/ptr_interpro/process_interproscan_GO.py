# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 12:24:32 2018

@author: weixiong
From interproscan output (already filtered for Pfam only), selects gene IDs with GO annotations and
writes them to a tab-delimited file, for uploading GO annotations with evidence codes to conekt
"""

import csv

original_file = 'ptr_aa_edited.fa.tsv'
new_file =  'GO_interproscan_pred.tsv'

with open(original_file, newline='') as original, open(new_file, 'w', newline='') as new:
    reader = csv.reader(original, delimiter='\t')
    writer = csv.writer(new, delimiter='\t')
    for line in reader:
        gene_ID = line[0]
        for item in line:
            if item.startswith('GO:'):
                GO_terms = item.split('|')
                for each in GO_terms:
                    new_line = [gene_ID, each, 'IEA']
                    writer.writerow(new_line)
