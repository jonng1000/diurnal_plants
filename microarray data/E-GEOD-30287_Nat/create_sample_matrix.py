# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 17:57:05 2018

@author: jonathan
"""
import csv
'''
1) this script creates a gene expression matrix .txt file, com_file.
2) it uses sample info and gene expression .txt files downloaded from Array
Express
3) sum_file is the .txt file showing which sample file is linked to which time
point
4) com_file is the .txt file containing the final matrix where each gene has
a list of gene expression values sorted according to time point, from earliest
to latest
'''
sum_file = 'E-GEOD-30287_summary.txt'
com_file = sum_file.split('_')[0] + '_combined.txt'

with open(sum_file, newline='') as file:
    filereader = csv.reader(file, delimiter='\t')
    next(filereader)
    sample_tp = [row for row in filereader]
    timepoints = [row[1] for row in sample_tp]  # used when writing file below

master_dict = {}
for row in sample_tp:
    sample_file = row[0] + '_sample_table.txt'
    with open(sample_file, newline='') as sfile:
        sfilereader = csv.reader(sfile, delimiter='\t')
        next(sfilereader)
        for row in sfilereader:
            if row[0] not in master_dict:
                # try except blocks to convert empty values (str) to 'NA'
                try:
                    master_dict[row[0]] = [float(row[1])]
                except ValueError:
                    if row[1] == '':
                        row[1] = 'NA'
                    master_dict[row[0]] = [row[1]]
            else:
                try:
                    master_dict[row[0]].append(float(row[1]))
                except ValueError:
                    if row[1] == '':
                        row[1] = 'NA'
                    master_dict[row[0]].append(row[1])

with open(com_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['Probe_ID'] + timepoints)
    for k, v in master_dict.items():
        writer.writerow([k] + v)

# this code block is to check to see that all rows in the matrix have the same
# length
# find the length of one key, then check that all keys have the
# same length
# if False is not printed, this implies that there are no missing values
one_gene = tuple(master_dict.keys())[0]
one_gene_values = len(master_dict[one_gene])
status = True
for k, v in master_dict.items():
    if len(v) == one_gene_values :
        pass
    else:
        print(k)
        print(False)
        status = False

# prints out success message if everything is fine
if status:
    print('Script ran without errors')
else:
    print('Please check, error encountered')