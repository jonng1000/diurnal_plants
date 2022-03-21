# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 11:46:18 2018

@author: weixiong
"""

from Bio import SeqIO

file1 = "Solanum_lycopersicum.SL2.50.cds.all.fa"
file2 = "ITAG3.2_CDS.fasta"

with open(file1) as f1, open(file2) as f2:
    records1 = SeqIO.parse(file1, 'fasta')
    records2 = SeqIO.parse(file2, 'fasta')
    
    dict1 = {}
    for record1 in records1:
        new_id1 = record1.id.split('.')[0]
        if new_id1 in dict1:
            raise ValueError('key already exist')
        else:
            dict1[new_id1] = record1.seq
        
    dict2 = {}
    for record2 in records2:
        new_id2 = record2.id.split('.')[0]
        if new_id2 in dict2:
            raise ValueError('key already exist')
        else:
            dict2[new_id2] = record2.seq
    
    true_dict1_dict2 = 0
    same_dict1_dict2 = 0
    not_same_dict1_dict2 = 0
    false_dict1_dict2 = 0        
    for key1 in dict1:
        if key1 in dict2:
            if dict1[key1] == dict2[key1]:
                same_dict1_dict2 += 1
            else:
                not_same_dict1_dict2 += 1
            true_dict1_dict2 += 1
        else:
            false_dict1_dict2 += 1
    total__dict1_dict2 = true_dict1_dict2 + false_dict1_dict2
    
    true_dict2_dict1 = 0
    same_dict2_dict1 = 0
    not_same_dict2_dict1 = 0
    false_dict2_dict1 = 0        
    for key2 in dict2:
        if key2 in dict1:
            if dict2[key2] == dict1[key2]:
                same_dict2_dict1 += 1
            else:
                not_same_dict2_dict1 += 1           
            true_dict2_dict1 += 1
        else:
            false_dict2_dict1 += 1
    total__dict2_dict1 = true_dict2_dict1 + false_dict2_dict1
    
    

