# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 20:16:55 2018

@author: weixiong

Compares the set of genes between two files and sees their intersection
"""

from Bio import Entrez, SeqIO
import pandas as pd

corrected_file = "ath_conekt.fa"
orig_gene_ex = 'JTK.Ath_Mat.txt'

record_dict = SeqIO.to_dict(SeqIO.parse(corrected_file, 'fasta'))
# length_keys is the number of genes from ensemble
length_keys = len(record_dict.keys()) # 27585
keys_set = set(record_dict.keys())
len(keys_set) == length_keys # True => no duplicates

gene_ex_df = pd.read_csv(orig_gene_ex, delimiter="\t", index_col=0)
gene_ex_index = gene_ex_df.index.values
# length_index is the number of genes from CF
length_index = len(gene_ex_index) # 22810
index_set = set(gene_ex_index)
len(index_set) == length_index # True => no duplicates
common = keys_set.intersection(index_set)
len_common = len(common) # 20977 genes identical between CF and Ensembl
percent_TAIR = len(common)/len(keys_set) * 100 # about 76%, about 6000 of CF's genes not in cds file

CF_not_cds = len(index_set - keys_set) # 1833 ensemble genes not in CF
cds_not_CF = len(keys_set - index_set) # 6608 ensemble genes not in CF

'''
# Ignore this code block, partially done as I wanted to use NCBI's functions
# to find gene synonyms for CF's genes which didn't appear in Ensembl, but MM
# said no need since only about 1000/22000 CF genes are not in Esenmbl

Entrez.email = "weixiong001@e.ntu.edu.sg"  

handle = Entrez.esearch(db="gene", term="At2g30510")
record = Entrez.read(handle)

if len(record["IdList"]) != 1:
    raise ValueError('More than one search term')
print(record)

handle = Entrez.efetch(db="gene", id="817601", retmode = 'text') #817601 works - gives one correct hit #At2g30510 doesnt work - gives multiple hits
print(handle.read())

record = SeqIO.read(handle, "fasta")
record = handler.read(handle)
handle.close()
print(record)
'''

