# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 12:05:28 2018

@author: weixiong

Takes interproscan output, and selects all rows from the Pfam database, and writes them
to a tab-delimited file, for uploading interproscan results to conekt
"""

import csv

original_file = 'ITAG3.2_proteins_interproscan.tsv'
new_file =  'interproscan_conekt.tsv'

with open(original_file, newline='') as original, open(new_file, 'w', newline='') as new:
    reader = csv.reader(original, delimiter='\t')
    writer = csv.writer(new, delimiter='\t')
    for line in reader:
        if line[3] == 'Pfam':
            gene_name = line[0].split('.')[0]
            line[0] = gene_name
            writer.writerow(line)
