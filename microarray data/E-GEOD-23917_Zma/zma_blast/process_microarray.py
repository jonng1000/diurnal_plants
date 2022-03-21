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
from Bio.Alphabet import IUPAC

pre_array = 'A-GEOD-10837_pre.txt'
new_array = 'zma_array_probes.txt'
fasta_blast = 'zma_array.fa'

with open(pre_array, newline='') as pre_a, open(new_array, 'w', newline='') as new_a:
    array_reader = csv.reader(pre_a, delimiter='\t')
    array_writer = csv.writer(new_a, delimiter='\t')
    for row in array_reader:
        if row[2]:
            array_writer.writerow([row[0], row[2]])
            
with open(new_array, newline='') as new_a, open(fasta_blast, 'w', newline='') as fb:
    array_reader = csv.reader(new_a, delimiter='\t')
    for row in array_reader:
        record = SeqRecord(Seq(row[1]), id=row[0], description=row[0])
        SeqIO.write(record, fb, 'fasta')
