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

original_file = "Arabidopsis_thaliana.TAIR10.pep.all.fa"
corrected_file = "ath_aa_conekt.fa"

with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = record.id.split('.')[0]
        isoform = record.id.split('.')[1]
        if isoform == '1':
            record.id = new_id
            # below removes the gene description, run either this or the line
            # further below to preserve description
            record.description = new_id 
            # below preserves the gene description
            #record.description = new_id + ' ' + record.description.split(' ', 1)[1]
            SeqIO.write(record, corrected, 'fasta')
        else:
            continue