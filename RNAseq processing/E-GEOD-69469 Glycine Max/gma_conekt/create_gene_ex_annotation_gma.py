# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 20:37:32 2018

@author: weixiong

Takes the gene expression matrix, original_file and and modifies it to make it
suitable for uploading into Conekt, by removing the decimals from the gene IDs
This produces the mod_file matrix

Also takes original_file and extracts out time points, and writes this plus
the time point's description, to an annotation file, new_file, which Conekt needs

However, duplicate descriptions in the annotation file will be treated as
replicates by Conekt, hence it needs to be manually checked to make sure that
this makes sense
"""

import csv

original_file = 'E-GEOD-69469_combined_PCA_edit.txt'
new_file =  'gene_ex_annotation.txt'
mod_file = 'gma_modified.txt'

with open(original_file, newline='') as original, open(new_file, 'w', newline='') as new, open(mod_file, 'w', newline='') as new2:
    reader = csv.reader(original, delimiter='\t')
    writer = csv.writer(new, delimiter='\t')
    row1 = next(reader)
    writer.writerow(['SampleID', 'Condition'])
    # creates the annotation file, new_file
    for time_pt in row1[1:]:
        without_dec = time_pt.split('.')[0]
        if without_dec[0] == 'L':
            writer.writerow([time_pt, without_dec[1:] + ' h after light'])
        elif without_dec[0] == 'D':
            writer.writerow([time_pt, without_dec[1:] + ' h after dark'])
        else:
            raise ValueError('Time point is neither L or D here')
    # creates the modified gene expression matrix, mod_file        
    writer2 = csv.writer(new2, delimiter='\t')
    writer2.writerow(row1)
    for line in reader:
        geneID = line[0].split('.')[0] + '.' + line[0].split('.')[1]
        mod_line = [geneID] + line[1:]
        writer2.writerow(mod_line)
