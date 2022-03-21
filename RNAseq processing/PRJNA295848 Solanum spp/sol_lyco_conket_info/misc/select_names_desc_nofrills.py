# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:54:29 2018

@author: weixiong

Creates a tab-delimited file with gene IDs (no decimals) in one column, and
gene description in another column
Does not use the Biopython module
"""
import csv

original_file = "ITAG3.2_CDS.fasta"
new_file = "name_description.txt" 
        
with open(original_file) as original, open(new_file, 'w', newline='') as new:
    for line in original:
        if line.startswith('>'):
            line = line[1:].rstrip('\n')
            raw_name = line.split(' ', 1)[0]
            name = raw_name.split('.')[0]
            try:
                desc = line.split(' ', 1)[1]
            except IndexError:
                desc = 'Function unknown'
        else:
            continue

        csvwriter = csv.writer(new, delimiter='\t')
        csvwriter.writerow([name, desc])