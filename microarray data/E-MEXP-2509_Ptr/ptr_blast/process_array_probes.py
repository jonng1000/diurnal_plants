# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 12:42:05 2019

@author: weixiong

Process ptr target genes and select their correct gene IDs
"""

import re
from Bio import SeqIO

original_file = "Poplar.target"
corrected_file = "ptr_array_probes.fa"

with open(original_file) as orig, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        new_id = re.split(':|;', record.description)[2]
        record.id = new_id
        record.description = new_id
        SeqIO.write(record, corrected, 'fasta')