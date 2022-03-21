# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:00:54 2018

@author: weixiong
"""

###############################################################################
# Creates a tab-delimited file which has both the transcript ID (old ID) and
# gene ID (new ID) in one line, corresponding to each gene
###############################################################################
from Bio import SeqIO
import csv
'''
original_file = "Nicotiana_attenuata.NIATTr2.cds.all.fa"
map_file = "new_old_IDs.txt"

with open(original_file) as orig, open(map_file, 'w', newline='') as mapf:
    records = SeqIO.parse(original_file, 'fasta')
    
    for record in records:
        old_id = record.id
        new_id = record.description.split('gene:')[1].split(' ')[0]        
        csvwriter = csv.writer(mapf, delimiter='\t')
        csvwriter.writerow([old_id, new_id])
'''
###############################################################################
# Modify original name_description.txt with new IDs, replaces transcript ID
# with gene ID, creates name_description_edited.txt
###############################################################################
'''
original_file = "name_description.txt"
edited_file = "name_description_edited.txt"
map_file = "new_old_IDs.txt"

with open(original_file) as orig, open(map_file, newline='') as mapf, \
     open(edited_file, 'w', newline='') as edit:
    reader1 = csv.reader(mapf, delimiter='\t')
    reader2 = csv.reader(orig, delimiter='\t')
    writer1 = csv.writer(edit, delimiter='\t')
    
    map_dict = {}
    for line in reader1:
        map_dict[line[1]] = line[0]
        
    for line in reader2:
        writer1.writerow([map_dict[line[0]], line[1]])
'''
###############################################################################
# Modify original name_description.txt with new IDs, replaces transcript ID
# with gene ID, creates name_description_edited.txt
###############################################################################
original_file = "nat_aa_conekt.fa.tsv"
edited_file = "nat_aa_conekt.fa_edited.tsv"
map_file = "new_old_IDs.txt"

with open(original_file) as orig, open(map_file, newline='') as mapf, \
     open(edited_file, 'w', newline='') as edit:
    reader1 = csv.reader(mapf, delimiter='\t')
    reader2 = csv.reader(orig, delimiter='\t')
    writer1 = csv.writer(edit, delimiter='\t')
    
    map_dict = {}
    for line in reader1:
        map_dict[line[1]] = line[0]
        
    for line in reader2:
        new_line = [map_dict[line[0]]] + line[1:]
        writer1.writerow(new_line)
        