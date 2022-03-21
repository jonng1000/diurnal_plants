# -*- coding: utf-8 -*-
"""
Spyder Editor
"""
'''
IMPT: THIS SCRIPT IS ONLY MEANT FOR DRA ACESSION NUMBERS 
Need to download CDS file (transciptome) first, as kallisto requires it.
'''   
import csv, os

sra_inputs = 'DRA006158_runs.txt'

# first, second, third and sra_accession are the variable names that when
# combined together, will form the complete donwload url
# need to write below string as r'...' as r indicates that it is a raw string
# so that escape characters will be ignored
first = 'D:/cli/bin/ascp.exe -i D:/cli/etc/asperaweb_id_dsa.openssh -QT -l 300m -P 33001 anonftp@ascp.ddbj.nig.ac.jp:/ddbj_database/dra/sra/ByExp/sra/DRX/DRX094/'

#D:/cli/bin/ascp.exe -i D:/cli/etc/asperaweb_id_dsa.openssh -QT -l 300m -P 33001 anonftp@ascp.ddbj.nig.ac.jp:/ddbj_database/dra/sra/ByExp/sra/DRX/DRX094/DRX094898/DRR101478/DRR101478.sra .
# name of index file output, produced by kallisto index
index_file = './brach_dist_transcripts.idx'
# name of cds file input, used by kallisto to produce .idx file
cds_file = './Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa.gz'

# produces index file, .idx from kallisto
#os.system('kallisto index -i ' + index_file + ' ' + cds_file)

with open(sra_inputs, newline='') as file:
    filereader = csv.reader(file, delimiter='\t')
    next(filereader)
    for row in filereader:
        second = row[0]
        sra_accession = row[1]
        dl_url = first + second + '/' + sra_accession + '/' + sra_accession + '.sra' + ' .'
        os.system(dl_url)
