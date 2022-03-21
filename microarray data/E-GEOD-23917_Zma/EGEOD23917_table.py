# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 17:57:05 2018

@author: jonathan
"""
import csv
'''
1) use sample info file from Array Express to create a file only containing
sample ID and timepoints
2) attach individual sample files according to timepoint
'''

with open('E-GEOD-23917_summary.txt', newline='') as file:
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
                # if the below code works, means that all values are numerical
                master_dict[row[0]] = [float(row[1])]
            else:
                master_dict[row[0]].append(float(row[1]))

# the length of one key is 18 values, hence below checks that all keys have the
# same length, implying that there are no missing values
for k, v in master_dict.items():
    if len(v) == 18:
        pass
    else:
        print(k)
        print(False)
        
with open('E-GEOD-23917_combined.txt', 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['Reporter Identifier'] + timepoints)
    for k, v in master_dict.items():
        writer.writerow([k] + v)
