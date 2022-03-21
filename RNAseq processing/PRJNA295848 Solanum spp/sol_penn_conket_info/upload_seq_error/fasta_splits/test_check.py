# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 17:35:58 2018

@author: weixiong

this script tests whether all the split files has the same number of genes as the
original, spe_conekt.fa

genes in spe_conekt.fa are counted here
"""

orig_counts = 0

for i in open('spe_conekt.fa','r').readlines():
    if '>' in i:
        orig_counts += 1
        
# conclusion: total number of genes in original file is 44965