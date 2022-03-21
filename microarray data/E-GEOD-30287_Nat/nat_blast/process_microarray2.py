# -*- coding: utf-8 -*-
"""
Spyder Editor

Takes A-GEOD-10837_pre.txt and produces zma_array.fa,
which are the probe sequences in fasta format for blast.
Produces zma_array_probes.txt as an intermediate file,
but it can be removed.
"""

import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

pre_array = 'A-GEOD-13527_pre.txt'
fasta_blast = 'nat_array.fa'

with open(pre_array, newline='') as pre_a, open(fasta_blast, 'w', newline='') as fb:
    array_reader = csv.reader(pre_a, delimiter='\t')
    for row in array_reader:
        if row[2]:
            record = SeqRecord(Seq(row[2]), id=row[0], description=row[0])
            SeqIO.write(record, fb, 'fasta')
