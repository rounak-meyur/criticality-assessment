# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:14:01 2020

@author: Rounak Meyur
Description: This program reads the dataset with bus data (name, id, kv and coordinates)
and .
"""


import sys,os
import networkx as nx
from pathlib import Path
from collections import namedtuple as nt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.sparse import csr_matrix
import numpy as np
import geopandas as gpd

workpath = os.getcwd()
datapath = str(Path(workpath).parent)+'/case/'
figpath = str(Path(workpath).parent)+'/figs/'

data_states = gpd.read_file(datapath+'states.shp')
state_polygon = list(data_states[data_states.STATE_ABBR == 
                                 'TX'].geometry.items())[0][1]

#%% Functions
def getbusdat(path,filename):
    """
    Extracts the bus data information from the txt file containing the information.

    Parameters
    ----------
    path : string
        path to the directory containing the file
    filename : string
        filename containing the bus data.

    Returns
    -------
    info : named tuple of type Busdata
        unimmutable record of bus data with fieldnames: cord, name and kv
        cord: geographical coordinates of the bus
        name: name of the bus location
        kv: base kV level of the bus

    """
    dict_cord = {}
    dict_name = {}
    dict_kv = {}
    with open(path+filename,'r') as file:
        for temp in file.readlines():
            data = temp.strip('\n').split('\t')
            dict_cord[int(data[0])] = [float(x) for x in data[3:]]
            dict_name[int(data[0])] = data[1]
            dict_kv[int(data[0])] = data[2]
    # Create named tuple with bus data
    bus = nt("Busdata",field_names=['cord','name','kv'])
    info = bus(cord=dict_cord,name=dict_name,kv=dict_kv)
    return info


def createnetwork(path,filename,criticality_file):
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
    e_count = {}
    
    # Get criticality values
    with open(path+criticality_file,'r') as file:
        crit_vals = [float(x) for x in file.readlines()[0].split(' ')]
    
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
            criticality[e] = crit_vals[i+1]
    graph = nx.MultiGraph()
    graph.add_edges_from(list(rate.keys()))
    nx.set_edge_attributes(graph,rate,'rating')
    nx.set_edge_attributes(graph,criticality,'criticality')
    return graph


def getgendata(path,filename):
    """
    Extracts generator information from required csv file.

    Parameters
    ----------
     path : string
        path to the directory containing the file
    filename : string
        filename containing the generator data.

    Returns
    -------
    dict_gen : dictionary object
        dictionary with keys as generator ID and values as bus ID and capacity of
        generation.
    condenser: list object
        list of synchronus condenser buses

    """
    dict_gen = {}
    condenser = []
    with open(path+filename,'r') as file:
        for i,temp in enumerate(file.readlines()[1:]):
            data = temp.strip('\n').split(',')
            if float(data[8])==0.0: condenser.append(int(data[0]))
            else: dict_gen[i+1]={'bus':int(data[0]),'cap':float(data[8])}
    return dict_gen,condenser


#%% Extract information
buses = getbusdat(datapath,'bus-dat.txt')
G = createnetwork(datapath,'branchdat.csv','criticality_1.5_2.txt')
GEN,synch_cond = getgendata(datapath,'gendat.csv')

genbus = [GEN[i]['bus'] for i in GEN]
gencap = [GEN[i]['cap'] for i in GEN]
nodelabel = {n:'G' if n in genbus else 'S' for n in list(G.nodes())}
for n in synch_cond: nodelabel[n] = 'C'
nx.set_node_attributes(G,nodelabel,'label')

# node attributes
nodelist = list(G.nodes())
colind = range(len(genbus))
rowind = [nodelist.index(b) for b in genbus]
data = [1 for _ in range(len(genbus))]
C = csr_matrix((data, (rowind, colind)), shape=(len(nodelist), len(genbus)))
cap = np.matmul(C.toarray(),np.array(gencap))

# edge attributes
edge_crit = nx.get_edge_attributes(G,'criticality')
edgelist = list(G.edges(keys=True))
top_critical = [e for e in edge_crit if edge_crit[e]<0.5]


#%% Plot the network
fig=plt.figure(figsize=(20,20))
ax=fig.add_subplot(111)

# Color and size nodes
nsize = []
ncolor = []
for i,n in enumerate(nodelist):
    if nodelabel[n]=='G':
        nsize.append(cap[i]/8.0)
        ncolor.append('lightsalmon')
    elif nodelabel[n] == 'C':
        nsize.append(20.0)
        ncolor.append('green')
    else:
        nsize.append(1.0)
        ncolor.append('blue')

ewidth = []
ecolor = []
for e in edgelist:
    if e in top_critical:
        ewidth.append(4.5)
        ecolor.append('crimson')
    else:
        ewidth.append(0.5)
        ecolor.append('black')


nx.draw_networkx(G,with_labels=False,ax=ax,pos=buses.cord,node_size=nsize,
                 node_color=ncolor,edgelist=edgelist,width=ewidth,style='dashed',
                 edge_color=ecolor,connectionstyle='arc3,rad=-3.0')
ax.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)

leglines = [Line2D([0], [0], color='black', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed'),
            Line2D([0], [0], color='crimson', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed'),
            Line2D([0], [0], color='white', markerfacecolor='lightsalmon', 
                   marker='o',markersize=10),
            Line2D([0], [0], color='white', markerfacecolor='blue', marker='o',
                   markersize=10),
            Line2D([0], [0], color='white', markerfacecolor='green', marker='o',
                   markersize=10)]

labels = ['transmission lines', 'top critical edges', 'generator buses', 
          'substation buses', 'synchronous condenser']

ax.legend(leglines,labels,loc='best',ncol=1,prop={'size': 15})
ax.set_title('Synthetic power grid of Texas with identified critical lines',
             fontsize=25)


for pol in state_polygon:
    x,y = pol.exterior.xy
    ax.plot(x,y,'g--')

figname = 'top-critical-1'
fig.savefig("{}{}.png".format(figpath,figname),bbox_inches='tight')










