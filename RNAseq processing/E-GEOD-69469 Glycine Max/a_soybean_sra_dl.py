# -*- coding: utf-8 -*-
"""
Spyder Editor
"""
'''
IMPT: THIS SCRIPT IS ONLY MEANT FOR SRA ACESSION NUMBERS WITH 7 DIGITS ONLY!!
example dl url:
D:\cli\bin\ascp.exe -QT -l 300m -P33001 -i D:\cli\etc\asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/
SRR245/004/SRR2452634 .

Need to download CDS file (transciptome) first, as kallisto requires it.
'''   
import csv, os

sra_inputs = 'soybean_two.txt'

# first, second, third and sra_accession are the variable names that when
# combined together, will form the complete donwload url
# need to write below string as r'...' as r indicates that it is a raw string
# so that escape characters will be ignored
first = r'D:\cli\bin\ascp.exe -QT -l 300m -P33001 -i D:\cli\etc\asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/'

# name of index file output, produced by kallisto index
index_file = r'.\soybean_transcripts.idx'
# name of cds file input, used by kallisto to produce .idx file
cds_file = r'.\Glycine_max.Glycine_max_v2.0.cds.all.fa.gz'

# produces index file, .idx from kallisto
os.system(r'.\kallisto index -i ' + index_file + ' ' + cds_file)

with open(sra_inputs, newline='') as file:
    filereader = csv.reader(file, delimiter='\t')
    for row in filereader:
        sra_accession = row[0]
        second = sra_accession[:6]
        last_digit = sra_accession[-1]
        third = '00' + last_digit
        dl_url = first + second + '/' + third + '/' + sra_accession + ' .'
        os.system(dl_url)
        # name of kallisto quant output folder
        out_file = sra_accession + '_kal_out'
        # list all files in dir then run below one by one
        os.system('./kallisto quant -i ' + index_file + ' -o ' + out_file +
                  ' --single -l 200 -s 20 ./SRR2048167/SRR2048167_1.fastq.gz