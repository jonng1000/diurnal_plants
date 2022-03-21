# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:00:54 2018

@author: weixiong
This script creates all the output for uploading files into conekt all at once.
Each code section marked by hashtags is created from individual smaller scripts,
I can run each section individually if necessary 
"""

###############################################################################
# From cds fasta file downloaded from the website, replaces headers with
# gene IDs (without decimal points), for uploading gene sequences to conekt 
###############################################################################
from Bio import SeqIO

original_file = "Solanum_lycopersicum.SL2.50.cds.all.fa"
corrected_file = "sly_conekt.fa"

with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = record.id.split('.')[0]
        record.id = new_id
        record.description = new_id
        SeqIO.write(record, corrected, 'fasta')
 
###############################################################################       
# Creates a tab-delimited file with gene IDs (no decimals) in one column, and
# gene description in another column
# Uses the Biopython module
# For uploading gene names and descriptions to conekt
###############################################################################
import csv

desc_file = "ITAG3.2_CDS.fasta"
gene_desc_file = "name_description.txt"

with open(desc_file) as orig_desc, open(gene_desc_file, 'w', newline='') as new_desc:
    records2 = SeqIO.parse(desc_file, 'fasta')
    for record2 in records2:
        name = record2.id.split('.')[0]
        try:
            desc = record2.description.split(' ', 1)[1]
        except IndexError:
            desc = 'Function unknown'

        csvwriter = csv.writer(new_desc, delimiter='\t')
        csvwriter.writerow([name, desc])

###############################################################################
# Takes interproscan output, and selects all rows from the Pfam database, and
# writes them to a tab-delimited file, for uploading interproscan results to 
# conekt
###############################################################################
interpro_output = 'ITAG3.2_proteins_interproscan.tsv'
pfam_file =  'interproscan_conekt.tsv'

with open(interpro_output, newline='') as orig_inter, open(pfam_file, 'w', newline='') as new_pfam:
    reader = csv.reader(orig_inter, delimiter='\t')
    writer = csv.writer(new_pfam, delimiter='\t')
    for line in reader:
        if line[3] == 'Pfam':
            gene_name = line[0].split('.')[0]
            line[0] = gene_name
            writer.writerow(line)

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
orig_gene_ex = 'PRJNA295848_sol_lyco_combined_PCA.txt'
new_file =  'gene_ex_annotation.txt'
mod_file = 'PRJNA295848_sol_lyco_modified.txt'

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

# Outputs file called gene_custers.txt
###############################################################################

input_f = 'JTK.sol_lyco.txt'
output_f = 'gene_custers.txt'

import numpy as np
import pandas as pd

gene_ex_df = pd.read_csv(input_f, delimiter="\t", index_col=0)
criteria = gene_ex_df['ADJ.P'] < 0.05
lower_than_p = gene_ex_df[criteria]

LAG_set = set(gene_ex_df['LAG'])
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