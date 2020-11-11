# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:14:01 2020

@author: Rounak Meyur
Description: This program reads the dataset with balance deviation values and aims
to plot the relative differences on the state map.
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
import pandas as pd

workpath = os.getcwd()
datapath = str(Path(workpath).parent)+'/case/'
resultpath = str(Path(workpath).parent)+'/results/'
figpath = str(Path(workpath).parent)+'/figs/'

data_states = gpd.read_file(datapath+'states.shp')
state_polygon = list(data_states[data_states.STATE_ABBR == 
                                 'TX'].geometry.items())[0][1]

#%% Functions
def get_neighbors(graph,edge,n=1):
    nodelist = [edge[0],edge[1]]
    new_nodes = [edge[0],edge[1]]
    for i in range(n):
        neighbors = []
        for node in new_nodes:
            neighbors.extend(list(nx.neighbors(graph,node)))
        neighbors = list(set(neighbors))
        new_nodes = [node for node in neighbors if node not in nodelist]
        nodelist = nodelist + new_nodes
    return nx.subgraph(graph,nodelist)



def getbusdat(path,filename,busfile):
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
        pd: scheduled load demand at bus

    """
    dict_cord = {}
    dict_name = {}
    dict_kv = {}
    dict_pd = {}
    with open(path+filename,'r') as file:
        for temp in file.readlines():
            data = temp.strip('\n').split('\t')
            dict_cord[int(data[0])] = [float(x) for x in data[3:]]
            dict_name[int(data[0])] = data[1]
            dict_kv[int(data[0])] = data[2]
    
    with open(path+busfile,'r') as file:
        for temp in file.readlines()[1:]:
            data = temp.strip('\n').split(',')
            dict_pd[int(data[0])] = float(data[2])
    
    # Create named tuple with bus data
    bus = nt("Busdata",field_names=['cord','name','kv','pd'])
    info = bus(cord=dict_cord,name=dict_name,kv=dict_kv,pd=dict_pd)
    return info


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
K = 1.5
buses = getbusdat(datapath,'bus-dat.txt','busdat.csv')
G = createnetwork(datapath,'branchdat.csv','criticality_base_'+str(K)+'_2.csv',
                  'balance_deviation.csv')
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

# Compute net power generation/consumption capacity
cap = np.matmul(C.toarray(),np.array(gencap))
load = np.array([buses.pd[n] for n in nodelist])
pinj = cap-load
pinj = {k:pinj[i] for i,k in enumerate(nodelist)}
nx.set_node_attributes(G,pinj,'power')

# edge attributes
edge_crit = nx.get_edge_attributes(G,'criticality')
edge_meanbaldev = nx.get_edge_attributes(G,'balanceMAD')
edge_medianbaldev = nx.get_edge_attributes(G,'balanceMedAD')
edgelist = list(G.edges(keys=True))
top_critical = [e for e in edge_crit if edge_crit[e]>0.75]
top_balance_dev = [e for e in edge_crit if edge_medianbaldev[e]>0.01]
rating = nx.get_edge_attributes(G,'rating')



#%% Plot the network
fig=plt.figure(figsize=(50,50))
ax=fig.add_subplot(111)

ewidth = []
ecolor = []
for e in edgelist:
    ewidth.append(int(edge_medianbaldev[e]))
    if e in top_balance_dev:
        # ewidth.append(10)
        ecolor.append('crimson')
    else:
        # ewidth.append(1)
        ecolor.append('black')

# Color/size by load and generation
nsize = []
ncolor = []
for n in nodelist:
    if pinj[n]>0.0:
        nsize.append(pinj[n]/2.0)
        ncolor.append('lightsalmon')
    elif pinj[n]<0.0:
        nsize.append(-pinj[n]/2.0)
        ncolor.append('royalblue')
    else:
        nsize.append(1.0)
        ncolor.append('limegreen')


nx.draw_networkx(G,with_labels=False,ax=ax,pos=buses.cord,node_size=nsize,
                 node_color=ncolor,edgelist=edgelist,width=ewidth,style='dashed',
                 edge_color=ecolor)
ax.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)

leglines = [Line2D([0], [0], color='black', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed'),
            Line2D([0], [0], color='crimson', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed'),
            Line2D([0], [0], color='white', markerfacecolor='lightsalmon', 
                   marker='o',markersize=25),
            Line2D([0], [0], color='white', markerfacecolor='royalblue', marker='o',
                   markersize=25),
            Line2D([0], [0], color='white', markerfacecolor='limegreen', marker='o',
                   markersize=25)]

labels = ['transmission lines', 'edges with high balance deviation', 
          'generator buses', 'load buses', 'zero power injection buses']

ax.legend(leglines,labels,loc='best',ncol=1,prop={'size': 25})
ax.set_title('Synthetic power grid of Texas with lines having high balance deviation',
             fontsize=25)


for pol in state_polygon:
    x,y = pol.exterior.xy
    ax.plot(x,y,'y--')

#%% 2nd figure

new_graph = nx.MultiGraph()
for e in top_balance_dev:
    sub_graph = get_neighbors(G,e,n=1)
    new_graph = nx.compose(new_graph,sub_graph)

nodelist = list(new_graph.nodes())
edgelist = list(new_graph.edges(keys=True))

fig=plt.figure(figsize=(50,50))
ax=fig.add_subplot(111)

ewidth = []
ecolor = []
for e in edgelist:
    if e in top_balance_dev:
        ewidth.append(500)
        ecolor.append('crimson')
    else:
        ewidth.append(1)
        ecolor.append('black')

# Color/size by load and generation
nsize = []
ncolor = []
for n in nodelist:
    if pinj[n]>0.0:
        nsize.append(pinj[n]/2.0)
        ncolor.append('lightsalmon')
    elif pinj[n]<0.0:
        nsize.append(-pinj[n]/2.0)
        ncolor.append('royalblue')
    else:
        nsize.append(1.0)
        ncolor.append('limegreen')



nx.draw_networkx(new_graph,with_labels=False,ax=ax,pos=buses.cord,node_size=0.0,
                 node_color=ncolor,edgelist=edgelist,width=ewidth,style='solid',
                 edge_color=ecolor,connectionstyle='arc3,rad=-3.0')
ax.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)

leglines = [Line2D([0], [0], color='black', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed'),
            Line2D([0], [0], color='crimson', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed'),
            Line2D([0], [0], color='white', markerfacecolor='lightsalmon', 
                   marker='o',markersize=25),
            Line2D([0], [0], color='white', markerfacecolor='royalblue', marker='o',
                   markersize=25),
            Line2D([0], [0], color='white', markerfacecolor='limegreen', marker='o',
                   markersize=25)]

labels = ['transmission lines', 'edges with high balance deviation', 
          'generator buses', 'load buses', 'zero power injection buses']

ax.legend(leglines,labels,loc='best',ncol=1,prop={'size': 25})
ax.set_title('Synthetic power grid of Texas with lines having high balance deviation',
             fontsize=25)


for pol in state_polygon:
    x,y = pol.exterior.xy
    ax.plot(x,y,'y--')


#%% Balance deviation causality
edgelist = list(G.edges(keys=True))
node_mismatch = [abs(sum([pinj[n]/100.0 for n in list(get_neighbors(G,e,n=0))])) \
                 for e in edgelist]
bal_dev = [edge_medianbaldev[e] for e in edgelist]
# node_bal_mismatch = [sum([pinj[n] for n in list(get_neighbors(G,e,n=0))]) \
#                      for e in top_balance_dev]
# bal_dev_2 = [edge_medianbaldev[e] for e in top_balance_dev]

num_comp = {}
for e in edgelist:
    net = G.copy()
    net.remove_edge(*e)
    num_comp[e] = nx.number_connected_components(net)


#%% Plot the causality
fig=plt.figure(figsize=(20,12))
ax=fig.add_subplot(111)
n1 = []; n2 = []
b1 = []; b2 = []

for i,e in enumerate(edgelist):
    if num_comp[e]==1:
        n1.append(node_mismatch[i])
        b1.append(bal_dev[i])
    else:
        n2.append(node_mismatch[i])
        b2.append(bal_dev[i])

ax.scatter(n1,b1,c='b',marker='+',s=50.0,label='Single connected component')
ax.scatter(n2,b2,c='r',marker='+',s=50.0,label='Isolated components/nodes')
ax.set_xlabel("Load-generation mismatch between disconnected nodes",fontsize=25)
ax.set_ylabel("Load-generation mismatch in the network",fontsize=25)
ax.legend(loc='best',prop={'size': 20})
fig.savefig("{}{}.png".format(figpath,'balance-dev-causality'),
            bbox_inches='tight')