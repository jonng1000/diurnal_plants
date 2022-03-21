# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os

# name of index file output, produced by kallisto index
index_file = 'sol_lyco_transcripts.idx'

for folder in os.listdir():
    if folder.startswith('SRR'):
        # name of kallisto quant output folder
        out_file = folder + '_kal_out'
        os.system('kallisto quant -i ' + index_file + ' -o ' + out_file +
                  ' --single -l 200 -s 20 ' + folder + '/' +
                  folder + '.fastq.gz')
