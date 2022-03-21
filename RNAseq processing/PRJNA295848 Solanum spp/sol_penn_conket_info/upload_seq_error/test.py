# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 17:35:58 2018

@author: weixiong

first, this script was ran by marek on spe_conekt.fa to change all multi
line sequences into a single line, then he modified it, such that its final
form is below, which replaces all 'N' in the sequence and checks that there
all sequence letters are in 'ATCGN'

produces spe_truncated.txt as output
"""

save = []

for i in open('spe_conekt.fa','r').readlines():
    if '>' in i:
        save.append('\n'+i)
    else:
        save.append(i.rstrip().replace('N','A'))
        if len(set(i.rstrip())-set('ATCGN'))>0:
            asdsa
        
v = open('spe_truncated.txt','w')
v.writelines(save)
v.close()