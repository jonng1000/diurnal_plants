# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:00:54 2018

@author: weixiong
"""

###############################################################################
# From cds fasta file downloaded from the website, replaces headers with
# gene IDs (without decimal points), for uploading gene sequences to conekt 
# modified from original create_conekt_files.py
# us this on both peptide and cds fasta files
###############################################################################
from Bio import SeqIO

#original_file = "Nicotiana_attenuata.NIATTr2.pep.all.fa"
#corrected_file = "nat_aa_conekt.fa"
original_file = "Nicotiana_attenuata.NIATTr2.cds.all.fa"
corrected_file = "nat_conekt.fa"

with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        
        # This is to create original nat_conekt.fa and nat_aa_conekt.fa files
        # which has correct IDs for Conekt, but Conekt doesn't seem to recognise
        # them. Remove if needed.
        #new_id = record.description.split('gene:')[1].split(' ')[0]
        #record.id = new_id
        #record.description = new_id
        
        # New code to use less correct IDs as Conekt accepts it
        # Comment out if not needed
        record.description = record.id
        
        SeqIO.write(record, corrected, 'fasta')
        
# Just to check if cds file has only unique genes
# Yes it has
record_dict = SeqIO.to_dict(SeqIO.parse(corrected_file, "fasta"))
test = [key for key in record_dict.keys()] # 33320 genes
test_set = set(test)
len(test_set) # 33320 genes
