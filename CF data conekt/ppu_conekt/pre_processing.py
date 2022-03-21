# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 22:03:33 2018

@author: weixiong

Takes the JTK output file and produces the gene_ex_annotation.txt and
gene_clusters.txt files. Also produces a modified.txt file which is supposed
to contain modified gene IDs, but for the input file, gene IDs do not need to be
modified
"""

import pandas as pd
import numpy as np
import csv

file = "JTK.Ppu_Mat_edit.txt"
'''
###############################################################################           
# Takes the gene expression matrix, original_file and modifies it to make it
# suitable for uploading into Conekt, by removing the decimals from the gene
# IDs. This produces the mod_file matrix

# Also takes original_file and extracts out time points, and writes this plus
# the time point's description, to an annotation file, new_file, which Conekt
# needs.

# However, duplicate descriptions in the annotation file will be treated as
# replicates by Conekt, hence it needs to be manually checked to make sure that
# this makes sense

# Modified from original create_conekt_files.py
###############################################################################
orig_gene_ex = file
new_file =  file.split('.')[1] + '_gene_ex_annotation.txt'
mod_file = file.split('.')[1] + '_modified.txt'

gene_ex_df = pd.read_csv(orig_gene_ex, delimiter="\t", index_col=0)
cleaned = gene_ex_df.drop(['BH.Q', 'ADJ.P', 'PER', 'AMP', 'LAG'], axis=1)

with open(new_file, 'w', newline='') as new:
    writer3 = csv.writer(new, delimiter='\t')
    writer3.writerow(['SampleID', 'Condition'])
    # creates the annotation file, new_file
    for time_pt in cleaned.columns.values:
        without_dec = time_pt.split('.')[0]
        if without_dec[0] == 'L':
            writer3.writerow([time_pt, without_dec[1:] + ' h after light'])
        elif without_dec[0] == 'D':
            writer3.writerow([time_pt, without_dec[1:] + ' h after dark'])
        else:
            raise ValueError('Time point is neither L or D here')

cleaned.to_csv(mod_file, sep='\t')
'''
###############################################################################
# Generate gene clusters, which are defined by genes with the same lag values.
# Sorts the cluster in ascending lag, where the first cluster is lag 0, second
# has a higher lag and so on.
# Also removes genes from CF JTK file which do not appear in the cds file.

# Outputs file called gene_clusters.txt
###############################################################################

input_f = file
output_f = file.split('.')[1] + '_gene_clusters.txt'
cds_file = "Porphyridium_CDS_UPDATED.fasta"

gene_ex_df = pd.read_csv(input_f, delimiter="\t", index_col=0)

# removes CF genes which does not appear in cds file
record_dict = SeqIO.to_dict(SeqIO.parse(cds_file, 'fasta'))
CF_genes = gene_ex_df.index.values
modified_genes = [x.split('|')[1] for x in CF_genes]
gene_ex_df.index = modified_genes
# not needed, but its just to see how many genes are there in common
keep_genes = [x for x in modified_genes if x in record_dict.keys()]
drop_genes = [x for x in modified_genes if x not in record_dict.keys()]
corrected = gene_ex_df.drop(drop_genes)

criteria = corrected['ADJ.P'] < 0.05
lower_than_p = corrected[criteria]

LAG_set = set(gene_ex_df['LAG'])
LAG_dict = {}

for period in LAG_set:
    selection = lower_than_p['LAG'] == period
    selected = lower_than_p[selection]
    selected_genes = selected.index
    LAG_dict[period] = selected_genes
    
LAG_dkeys = list(LAG_dict.keys())
LAG_dkeys.sort()

with open(output_f, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    for k in range(24): #LAG_dkeys # also put a blank list below
        if k in LAG_dict:
            writer.writerow(LAG_dict[k])
        else:
            writer.writerow(['abc'])
