# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:00:54 2018

@author: weixiong

From cds file, remove decimals as they indicate splicing isoforms. Only gene
IDs with a .1 at the end are copied, along with their fasta sequences into
corrected_file, which is the fasta file to be uploaded into Conekt. Gene IDs
with a .2 or above are ignored
"""

###############################################################################
# From cds fasta file downloaded from the website, replaces headers with
# gene IDs (without decimal points), for uploading gene sequences to conekt 
###############################################################################
from Bio import SeqIO

original_file = "Spenn-v2-cds-annot.fa"
corrected_file = "spe_conekt.fa"

with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    
    # this is just to see how many unique genes there are after splicing
    # isoforms are removed
    # records_dict = SeqIO.to_dict(records)
    # 48923 genes
    # all_keys = records_dict.keys()
    # 44965 after removing isoforms
    # unique_genes = {gene.split('.')[0] for gene in all_keys}  

    # if i count number of times genes are written to file using below code,
    # i get 44965, same as above   
    for record in records:
        new_id = record.id.split('.')[0]
        isoform = record.id.split('.')[1]
        if isoform == '1':
            record.id = new_id
            # below removes the gene description, run either this or the line
            # further below to preserve description
            record.description = new_id 
            # below preserves the gene description
            #new_id + ' ' + record.description.split(' ', 1)[1]
            SeqIO.write(record, corrected, 'fasta')
        else:
            continue