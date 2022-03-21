# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import csv
'''
Sample timepoint info for this dataset shows time points linked to expt accession
IDs, not run IDs, hence this script replaces ext accession IDs with run IDs
for the timepoint info
'''

timepts = 'timepts_expt_ID.txt'
runs_info = 'PRJNA263611_runs.txt'
sample_dict = {}

with open(runs_info, newline='') as f:
    freader = csv.reader(f, delimiter='\t')
    next(freader)
    for r in freader:
        if r[0] not in sample_dict:
            sample_dict[r[0]] = r[1]

with open(timepts, newline='') as file:
    filereader = csv.reader(file, delimiter='\t')
    header = next(filereader)
    table_body = [row for row in filereader]
    com_table = [header] + table_body

with open('draft_PRJNA263611_summary.txt', 'w', newline='') as file1:
    writer = csv.writer(file1, delimiter='\t')
    writer.writerow(com_table[0])
    for each_row in com_table[1:]:
        sample_IDs = [sample_dict[item] for item in each_row]
        writer.writerow(sample_IDs)
