# -*- coding: utf-8 -*-
"""
Spyder Editor
"""
'''
IMPT: THIS SCRIPT IS ONLY MEANT FOR SRA ACESSION NUMBERS WITH 6 OR 7 DIGITS ONLY!!
example dl url:
D:\cli\bin\ascp.exe -QT -l 300m -P33001 -i D:\cli\etc\asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/
SRR245/004/SRR2452634 .

Need to download CDS file (transciptome) first, as kallisto requires it.
It can be used for when seq data is either paired or single reads

Modified from pseu_men_sra_dl.py in the Windows desktop, so that it in can run in                               
the Linux workstation.
'''   
import csv, os, shutil

sra_inputs = 'PRJNA263611_runs.txt'

# first, second, third and sra_accession are the variable names that when
# combined together, will form the complete donwload url
# need to write below string as r'...' as r indicates that it is a raw string
# so that escape characters will be ignored
first = 'ascp -QT -l 300m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/'

# name of index file output, produced by kallisto index
index_file = 'pseu_men_transcripts.idx'
# name of cds file input, used by kallisto to produce .idx file
cds_file = 'cds.pme.csv.gz'

# produces index file, .idx from kallisto
#os.system('kallisto index -i ' + index_file + ' ' + cds_file)

with open(sra_inputs, newline='') as file:
    filereader = csv.reader(file, delimiter='\t')
    next(filereader)
    for row in filereader:
        sra_accession = row[1]
        single_paired = row[2]  # paired end or single end seq
        second = sra_accession[:6]
        last_digit = sra_accession[-1]
        third = '00' + last_digit
        # run accession no has 6 digits, but has length 9 due to 3 letters in
        # front
        if len(sra_accession) == 9:
            dl_url = first + second + '/' + sra_accession + ' .'
        # run accession no has 7 digits, but has length 9 due to 3 letters in
        # front
        elif len(sra_accession) == 10:
            dl_url = first + second + '/' + third + '/' + sra_accession + ' .'
        os.system(dl_url)
        # name of kallisto quant output folder
        out_file = sra_accession + '_kal_out'
        if single_paired == 'SINGLE':
            os.system('kallisto quant -i ' + index_file + ' -o ' + out_file +
                      ' --single -l 200 -s 20 ./' + sra_accession + '/' +
                      sra_accession + '.fastq.gz')
        elif single_paired == 'PAIRED':
            os.system('kallisto quant -i ' + index_file + ' -o ' + out_file +
                      ' --single -l 200 -s 20 ./' + sra_accession + '/' +
                      sra_accession + '_1.fastq.gz')
        #shutil.rmtree('./' + sra_accession)
