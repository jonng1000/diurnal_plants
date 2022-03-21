# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 17:35:58 2018

@author: weixiong

this script tests whether all the split files has the same number of genes as the
original, spe_conekt.fa

genes in the split spe_conekt.fa files are counted here
"""
import os

split_counts = 0
for f in os.listdir():
    if f.startswith('spe'):
        file =  open(f,'r')

        for i in file.readlines():
            if '>' in i:
                split_counts += 1
                
# conclusion: total number of genes in split files is 44965