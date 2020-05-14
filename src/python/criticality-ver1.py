# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:14:01 2020

@author: Rounak Meyur
Description: This program reads the dataset with bus coordinates and creates a txt file
with useful information and displayed in a systematic manner.
"""


import sys,os
import networkx as nx
from pathlib import Path

workpath = os.getcwd()
datapath = str(Path(workpath).parent)+'/case/'

D = {}
with open(datapath+'texas-busdat.txt','r') as file:
    lines = file.readlines()
    for temp in lines:
        data = temp.strip('\n').split('      ')
        name = data[0].split('"')[1]
        if data[0][-1]=='"':
            kv = '%0.2f'%(float(data[1].split('"')[0]))
        else:
            kv = '%0.2f'%(float(data[0].split('"')[2]))
        cord = [float(x) for x in data[-1].split('    ')[1].split('"')[0].strip(' ').split('  ')][::-1]
        D[int(data[0][:4])]={'cord':cord,'name':name,'kv':kv}


refined = '\n'.join(['\t'.join([str(d),D[d]['name'],D[d]['kv'],
                                str(D[d]['cord'][0]),str(D[d]['cord'][1])]) for d in D])

with open(datapath+'bus-dat.txt','w') as f:
    f.write(refined)


#%% Criticality data
with open(datapath+'criticality_1.5_2.txt','r') as file:
    lines = ' '.join([temp.strip('\n') for temp in file.readlines()])
    vals = lines.split(' ')

with open(datapath+'criticality_1.5_2.txt','w') as file:
    file.write(lines)