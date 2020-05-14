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

workpath = os.getcwd()
datapath = str(Path(workpath).parent)+'/case/'



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


def createnetwork(path,filename):
    """
    Creates a networkx multigraph representing the network to be analyzed.

    Parameters
    ----------
    path : string
        path to the directory containing the file
    filename : string
        filename containing the branch data.

    Returns
    -------
    graph : networkx undirected multigraph
        Networkx graph with edge attributes as the rating of the lines.

    """
    rate = {}
    e_count = {}
    with open(path+filename,'r') as file:
        for temp in file.readlines()[1:]:
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
    graph = nx.MultiGraph()
    graph.add_edges_from(list(rate.keys()))
    nx.set_edge_attributes(graph,rate,'rating')
    return graph



buses = getbusdat(datapath,'bus-dat.txt')
G = createnetwork(datapath,'branchdat.csv')


#%% Plot the points
# xpt = [d[0] for d in list(D.values())]
# ypt = [d[1] for d in list(D.values())]
# import matplotlib.pyplot as plt
# fig=plt.figure(figsize=(20,20))
# ax=fig.add_subplot(111)
# ax.scatter(xpt,ypt,s=10.0,c='r')