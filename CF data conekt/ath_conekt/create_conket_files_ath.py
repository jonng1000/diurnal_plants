# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:00:54 2018

@author: weixiong
Modified from create_conket_files.py, which creates all the output for
uploading files into conekt all at once. Each code section marked by hashtags
is created from individual smaller scripts, I can run each section
individually if necessary. Some sections have been modified from the original,
create_conket_files.py 
"""

###############################################################################
# From cds fasta file downloaded from the website, replaces headers with
# gene IDs (without decimal points), for uploading gene sequences to conekt 
###############################################################################
from Bio import SeqIO

# this section is deleted as another script was ran instead
 
###############################################################################       
# Creates a tab-delimited file with gene IDs (no decimals) in one column, and
# gene description in another column
# Uses the Biopython module
# For uploading gene names and descriptions to conekt
###############################################################################
import csv

# this section is deleted as another script was ran instead

###############################################################################
# Takes interproscan output, and selects all rows from the Pfam database, and
# writes them to a tab-delimited file, for uploading interproscan results to 
# conekt
###############################################################################

pfam_file =  'ath_aa_conekt_edited.fa.tsv'

# this section is deleted as interproscan output only had pfam results

###############################################################################
# From interproscan output (already filtered for Pfam only), selects gene IDs
# with GO annotations and writes them to a tab-delimited file, for uploading
# GO annotations with evidence codes to conekt
###############################################################################

GO_file =  'GO_interproscan_pred.tsv'

with open(pfam_file, newline='') as orig_GO, open(GO_file, 'w', newline='') as new_GO:
    reader2 = csv.reader(orig_GO, delimiter='\t')
    writer2 = csv.writer(new_GO, delimiter='\t')
    for line in reader2:
        gene_ID = line[0]
        for item in line:
            if item.startswith('GO:'):
                GO_terms = item.split('|')
                for each in GO_terms:
                    new_line = [gene_ID, each, 'IEA']
                    writer2.writerow(new_line)
exit()
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
###############################################################################
orig_gene_ex = 'PRJNA295848_sol_penn_combined_PCA.txt'
new_file =  'gene_ex_annotation.txt'
mod_file = 'PRJNA295848_sol_penn_modified.txt'

with open(orig_gene_ex, newline='') as original, open(new_file, 'w', newline='') as new, open(mod_file, 'w', newline='') as new2:
    reader3 = csv.reader(original, delimiter='\t')
    writer3 = csv.writer(new, delimiter='\t')
    row1 = next(reader3)
    writer3.writerow(['SampleID', 'Condition'])
    # creates the annotation file, new_file
    for time_pt in row1[1:]:
        without_dec = time_pt.split('.')[0]
        if without_dec[0] == 'L':
            writer3.writerow([time_pt, without_dec[1:] + ' h after light'])
        elif without_dec[0] == 'D':
            writer3.writerow([time_pt, without_dec[1:] + ' h after dark'])
        else:
            raise ValueError('Time point is neither L or D here')
    # creates the modified gene expression matrix, mod_file        
    writer4 = csv.writer(new2, delimiter='\t')
    writer4.writerow(row1)
    for line in reader3:
        geneID = line[0].split('.')[0]
        mod_line = [geneID] + line[1:]
        writer4.writerow(mod_line)

###############################################################################
# Generate gene clusters, which are defined by genes with the same lag values.
# Sorts the cluster in ascending lag, where the first cluster is lag 0, second
# has a higher lag and so on.

# Outputs file called gene_clusters.txt
###############################################################################

input_f = 'JTK.Ath_Mat.txt'
cds_file = 'ath_conekt.fa'
output_f = 'gene_clusters.txt'

import numpy as np
import pandas as pd

gene_ex_df = pd.read_csv(input_f, delimiter="\t", index_col=0)
record_dict = SeqIO.to_dict(SeqIO.parse(cds_file, 'fasta'))
# removes CF genes which do not appear in cds file
# not needed, but its just to see how many genes are there in common
keep_genes = [x for x in gene_ex_df.index.values if x in record_dict.keys()]
drop_genes = [x for x in gene_ex_df.index.values if x not in record_dict.keys()]
corrected = gene_ex_df.drop(drop_genes)

criteria = corrected['ADJ.P'] < 0.05
lower_than_p = corrected[criteria]

LAG_set = set(corrected['LAG'])
LAG_dict = {}

for period in LAG_set:
    selection =lower_than_p['LAG'] == period
    selected =lower_than_p[selection]
    selected_genes = [item.split('.')[0] for item in selected.index]
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