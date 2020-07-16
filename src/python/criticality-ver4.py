# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 22:16:28 2020

@author: meyu507
"""


import os
from pathlib import Path
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt


workpath = os.getcwd()
datapath = str(Path(workpath).parent)+'/case/'
resultpath = str(Path(workpath).parent)+'/results/'
figpath = str(Path(workpath).parent)+'/figs/'

cont = [85,87,464,3007]
fig = plt.figure(figsize=(20,20))

for i,c in enumerate(cont):
    ax = fig.add_subplot(4,1,i+1)
    filename = 'results_'+str(c)+'.csv'
    D = genfromtxt(resultpath+filename, delimiter=',',skip_header=1)[:,1:]
    D[c-1,:] = np.zeros(shape=(1,1000))
    S = np.where(np.sum(D,1)>100)
    ax.spy(D.T,markersize=2,color='red')
