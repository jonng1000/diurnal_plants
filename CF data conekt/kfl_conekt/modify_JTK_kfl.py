# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:55:22 2019

@author: weixiong

This script removes gemes from CF's data which do not appear in the cds file,
and writes a new tab-delimited file with all the genes that are left. Modified
from modify_JTK.py.
"""
import pandas as pd
from Bio import SeqIO

orig_gene_ex = 'JTK.Kfl_Mat_edit.txt'
cds_file = "160614_klebsormidium_v1.1_CDS.fasta"
correct_file = 'Kfl_Mat_edit_modified300119.txt'

gene_ex_df = pd.read_csv(orig_gene_ex, delimiter="\t", index_col=0)
record_dict = SeqIO.to_dict(SeqIO.parse(cds_file, 'fasta'))

edited_df = gene_ex_df.drop(['BH.Q', 'ADJ.P', 'PER', 'LAG', 'AMP'], axis=1)
CF_genes = gene_ex_df.index.values

# not needed, but its just to see how many genes are there in common
keep_genes = [x for x in CF_genes if x in record_dict.keys()]
drop_genes = [x for x in CF_genes if x not in record_dict.keys()]
corrected = edited_df.drop(drop_genes)

corrected.to_csv(correct_file, sep='\t')
