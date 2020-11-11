# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 22:16:28 2020

@author: meyu507
"""


import os,sys
from pathlib import Path
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx


workpath = os.getcwd()
datapath = str(Path(workpath).parent)+'/case/'
resultpath = str(Path(workpath).parent)+'/results/'
figpath = str(Path(workpath).parent)+'/figs/'


def get_edgelist(path,filename):
    e_count = {}
    edgelist = []
    with open(path+filename,'r') as file:
        for i,temp in enumerate(file.readlines()[1:]):
            data = temp.strip('\n').split(',')
            edge = (int(data[0]),int(data[1]))
            if edge not in e_count and tuple(list(edge)[::-1]) not in e_count:
                e_count[edge] = 1
                e = (int(data[0]),int(data[1]),str(e_count[edge]))
            elif edge in e_count: 
                e_count[edge] += 1
                e = (int(data[0]),int(data[1]),str(e_count[edge]))
            else: 
                e_count[tuple(list(edge)[::-1])] += 1
                e = (int(data[0]),int(data[1]),str(e_count[tuple(list(edge)[::-1])]))
            edgelist.append(e)
    return edgelist

def createnetwork(path,filename,criticality_file,balancedev_file):
    """
    Creates a networkx multigraph representing the network to be analyzed.

    Parameters
    ----------
    path : string
        path to the directory containing the file
    filename : string
        filename containing the branch data.
    criticality_file : string
        filename containing the criticality data.

    Returns
    -------
    graph : networkx undirected multigraph
        Networkx graph with edge attributes as the rating of the lines.

    """
    
    rate = {}
    criticality = {}
    mean_baldev = {}
    median_baldev = {}
    e_count = {}
    
    # Get criticality values
    with open(path+criticality_file,'r') as file:
        line_file = file.readlines()[2:]
        crit_vals = [float(x.strip('\n').split(',')[1]) for x in line_file]
    
    # Get criticality values
    with open(path+balancedev_file,'r') as file:
        line_file = file.readlines()[1:]
        mean_baldev_vals = [float(x.strip('\n').split(',')[0]) for x in line_file]
        median_baldev_vals = [float(x.strip('\n').split(',')[1]) for x in line_file]
    
    # Get edgelist and ratings
    with open(path+filename,'r') as file:
        for i,temp in enumerate(file.readlines()[1:]):
            data = temp.strip('\n').split(',')
            edge = (int(data[0]),int(data[1]))
            if edge not in e_count and tuple(list(edge)[::-1]) not in e_count:
                e_count[edge] = 1
                e = (int(data[0]),int(data[1]),str(e_count[edge]))
            else:
                if edge in e_count: 
                    e_count[edge] += 1
                    e = (int(data[0]),int(data[1]),str(e_count[edge]))
                else: 
                    e_count[tuple(list(edge)[::-1])] += 1
                    e = (int(data[0]),int(data[1]),str(e_count[tuple(list(edge)[::-1])]))
            rate[e] = float(data[5])
            criticality[e] = crit_vals[i]
            mean_baldev[e] = mean_baldev_vals[i]
            median_baldev[e] = median_baldev_vals[i]
    graph = nx.MultiGraph()
    graph.add_edges_from(list(rate.keys()))
    nx.set_edge_attributes(graph,rate,'rating')
    nx.set_edge_attributes(graph,criticality,'criticality')
    nx.set_edge_attributes(graph,mean_baldev,'balanceMAD')
    nx.set_edge_attributes(graph,median_baldev,'balanceMedAD')
    return graph



K = 1.5
G = createnetwork(datapath,'branchdat.csv','criticality_base_'+str(K)+'_2.csv',
                  'balance_deviation.csv')
edge_crit = nx.get_edge_attributes(G,'criticality')
all_edges = get_edgelist(datapath,'branchdat.csv')
top_critical = [e for e in edge_crit if edge_crit[e]>0.75]
cont = sorted([all_edges.index(edge) if edge in all_edges \
        else all_edges.index((edge[1],edge[0],edge[2])) \
            for edge in top_critical])


figfile = 'overloaded-'+str(K)
overloaded = {}

# cont = [85,87,464,3007]
# cont = [50,85,87,93,114,116,400,432,464,1773,1774,2058,2101,3007]
# fig = plt.figure(figsize=(20,20))

#%% Various thresholds
thresh = 500
for i,c in enumerate(cont):
    # ax = fig.add_subplot(2,2,i+1)
    filename = 'results_'+str(K)+'_'+str(c+1)+'.csv'
    D = genfromtxt(resultpath+filename, delimiter=',',skip_header=1)
    D[c,:] = np.zeros(shape=(1,1000))
    ind = np.where(np.sum(D,1)>thresh)[0].tolist()
    overloaded[c+1] = [j+1 for j in ind]
    # ax.spy(D.T,markersize=2,color='red')
    # ax.tick_params(left=False,top=False,labelleft=False,labeltop=False)
    # ax.set_ylabel("Scenarios")
    # ax.set_xlabel("Lines")
    # fig.savefig("{}{}.png".format(figpath,figfile),bbox_inches='tight')

print(overloaded)


