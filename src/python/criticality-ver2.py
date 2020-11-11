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
import pandas as pd

workpath = os.getcwd()
datapath = str(Path(workpath).parent)+'/case/'
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
        line_file = file.readlines()[2:]
        crit_vals = [float(x.strip('\n').split(',')[1]) for x in line_file]
    
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


#%% Extract information
K = 1.5
buses = getbusdat(datapath,'bus-dat.txt','busdat.csv')
G = createnetwork(datapath,'branchdat.csv','criticality_base_'+str(K)+'_2.csv')
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
# edgelist = list(G.edges(keys=True))
all_edges = get_edgelist(datapath,'branchdat.csv')
top_critical = [e for e in edge_crit if edge_crit[e]>0.75]
cont = sorted([all_edges.index(edge) if edge in all_edges \
        else all_edges.index((edge[1],edge[0],edge[2])) \
            for edge in top_critical])
rating = nx.get_edge_attributes(G,'rating')

#%% Data obtained from R code
## top LODF lines
top_lodf = {84:[86,115,33],86:[84,34,65],
            463:[460,2046,2047],3006:[3154,3167,3004]}
top_olf = {84:[34,33,57],86:[34,33,57],
            463:[460,34,2135],3006:[3047,3073,34]}
top_scenario = {84:[33,34],86:[33,34,57],463:[34,460],3006:[3046,3047,3073]}

i = 3
z = 1.2
a = 1
m=5
cont_id = cont[i]
critical = all_edges[cont_id]
t_lodf = [all_edges[j] for j in top_lodf[cont_id]]
t_olf = [all_edges[j] for j in top_olf[cont_id]]
t_sce = [all_edges[j] for j in top_scenario[cont_id]]

#%% Plot the network
fig=plt.figure(figsize=(50,50))
ax=fig.add_subplot(111)

ewidth = []
ecolor = []
for e in all_edges:
    alt_e = (e[1],e[0],e[2])
    if (e == critical) or (alt_e == critical):
        # ewidth.append(rating[e]/1000.0)
        ewidth.append(10)
        ecolor.append('crimson')
    # elif (e in t_olf) or (alt_e in t_olf):
    #     ewidth.append(10)
    #     ecolor.append('orchid')
    # elif (e in t_lodf) or (alt_e in t_lodf):
    #     ewidth.append(10)
    #     ecolor.append('darkgreen')
    else:
        # ewidth.append(rating[e]/1000.0)
        ewidth.append(0.1)
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
                 node_color=ncolor,edgelist=all_edges,width=ewidth,style='dashed',
                 edge_color=ecolor,connectionstyle='arc3,rad=-3.0')
ax.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)

leglines = [Line2D([0], [0], color='crimson', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed',linewidth=10),
            Line2D([0], [0], color='white', markerfacecolor='lightsalmon', 
                   marker='o',markersize=25),
            Line2D([0], [0], color='white', markerfacecolor='royalblue', marker='o',
                   markersize=25),
            Line2D([0], [0], color='black', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed',linewidth=10),
            Line2D([0], [0], color='green', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed',linewidth=10),
            Line2D([0], [0], color='orchid', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed',linewidth=10),
            Line2D([0], [0], color='blue', markerfacecolor='white', marker='*',
                   markersize=0,linestyle='dashed',linewidth=10)]

labels = ['line outage contingency', 'generator buses', 'load buses',
          'transmission lines','lines with high LODF',
          'lines with high overload factor','overloaded lines in random scenarios']

ax.legend(leglines,labels,loc='best',ncol=1,prop={'size': 45})
ax.set_title('Synthetic power grid of Texas with identified critical lines',
             fontsize=50)


for pol in state_polygon:
    x,y = pol.exterior.xy
    ax.plot(x,y,'y--')




#%% Examine individual contingency

def draw_subnetwork_top(ax,graph,edge_interest,n=4):
    """
    Creates the sub network from the original network which just includes the
    neighbors of the interesting edge. The figure displays the rating of the
    lines in the neighborhood of the critical line.

    Parameters
    ----------
    ax : matplotlib axes object
        The axes on which the networkx network would be plotted.
    graph: networkx MultiGraph object
        The networkx graph representing the power grid network.
    n : integer, optional
        Number of hops to consider from the interesting edge.. The default is 4.

    Returns
    -------
    ax : matplotlib axes object
        The axes on which the networkx network would be plotted.

    """
    rate = nx.get_edge_attributes(graph,'rating')
    p_inj = nx.get_node_attributes(graph,'power')
    interest = get_neighbors(graph,edge_interest,n)
    xpts = [buses.cord[n][0] for n in list(interest.nodes())]
    ypts = [buses.cord[n][1] for n in list(interest.nodes())]
    
    # Color by rating
    edgelist_inset = list(interest.edges(keys=True))
    nodelist_inset = list(interest.nodes())
    ewidth = []
    ecolor = []
    for e in edgelist_inset:
        alt_e = (e[1],e[0],e[2])
        if e == critical or alt_e == critical:
            ewidth.append(10)
            # if e in rate: ewidth.append(rate[e]/200.0)
            # else: ewidth.append(rate[alt_e]/200.0)
            ecolor.append('crimson')
        elif (e in t_olf) or (alt_e in t_olf):
            ewidth.append(10)
            ecolor.append('orchid')
        elif (e in t_lodf) or (alt_e in t_lodf):
            ewidth.append(10)
            ecolor.append('darkgreen')
        else:
            ewidth.append(3)
            # if e in rate: ewidth.append(rate[e]/200.0)
            # else: ewidth.append(rate[alt_e]/200.0)
            ecolor.append('black')
    
    # Color/size by load and generation
    nsize = []
    ncolor = []
    for n in nodelist_inset:
        if p_inj[n]>0.0:
            nsize.append(p_inj[n]/1.0)
            ncolor.append('lightsalmon')
        elif p_inj[n]<0.0:
            nsize.append(-p_inj[n]/1.0)
            ncolor.append('royalblue')
        else:
            nsize.append(2.0)
            ncolor.append('limegreen')
    
    # elabel={(r[0],r[1]):rate[r] for r in rate \
    #         if r in edgelist_inset or (r[1],r[0],r[2]) in edgelist_inset}
    nx.draw_networkx(interest,with_labels=False,ax=ax,nodelist=nodelist_inset,
                     pos=buses.cord,node_size=nsize,node_color=ncolor,
                     edgelist=edgelist_inset,width=ewidth,style='dashed',
                     edge_color=ecolor)
    # nx.draw_networkx_edge_labels(interest,ax=ax,pos=buses.cord,edge_labels=elabel,
    #                              font_weight='normal',font_size=25)
    
    
    ax.set_xlim(min(xpts),max(xpts))
    ax.set_ylim(min(ypts),max(ypts))
    ax.tick_params(bottom=False,left=False,labelleft=False,labelbottom=False)
    return ax


def draw_subnetwork_bottom(ax,graph,edge_interest,n=4):
    """
    Creates the sub network from the original network which just includes the
    neighbors of the interesting edge. The figure displays the rating of the
    lines in the neighborhood of the critical line.

    Parameters
    ----------
    ax : matplotlib axes object
        The axes on which the networkx network would be plotted.
    graph: networkx MultiGraph object
        The networkx graph representing the power grid network.
    n : integer, optional
        Number of hops to consider from the interesting edge.. The default is 4.

    Returns
    -------
    ax : matplotlib axes object
        The axes on which the networkx network would be plotted.

    """
    rate = nx.get_edge_attributes(graph,'rating')
    p_inj = nx.get_node_attributes(graph,'power')
    interest = get_neighbors(graph,edge_interest,n)
    xpts = [buses.cord[n][0] for n in list(interest.nodes())]
    ypts = [buses.cord[n][1] for n in list(interest.nodes())]
    
    # Color by rating
    edgelist_inset = list(interest.edges(keys=True))
    nodelist_inset = list(interest.nodes())
    ewidth = []
    ecolor = []
    for e in edgelist_inset:
        alt_e = (e[1],e[0],e[2])
        if e == critical or alt_e == critical:
            ewidth.append(10)
            # if e in rate: ewidth.append(rate[e]/200.0)
            # else: ewidth.append(rate[alt_e]/200.0)
            ecolor.append('crimson')
        elif (e in t_sce) or (alt_e in t_sce):
            ewidth.append(10)
            ecolor.append('blue')
        else:
            ewidth.append(3)
            # if e in rate: ewidth.append(rate[e]/200.0)
            # else: ewidth.append(rate[alt_e]/200.0)
            ecolor.append('black')
    
    # Color/size by load and generation
    nsize = []
    ncolor = []
    for n in nodelist_inset:
        if p_inj[n]>0.0:
            nsize.append(p_inj[n]/1.0)
            ncolor.append('lightsalmon')
        elif p_inj[n]<0.0:
            nsize.append(-p_inj[n]/1.0)
            ncolor.append('royalblue')
        else:
            nsize.append(2.0)
            ncolor.append('limegreen')
    
    # elabel={(r[0],r[1]):rate[r] for r in rate \
    #         if r in edgelist_inset or (r[1],r[0],r[2]) in edgelist_inset}
    nx.draw_networkx(interest,with_labels=False,ax=ax,nodelist=nodelist_inset,
                     pos=buses.cord,node_size=nsize,node_color=ncolor,
                     edgelist=edgelist_inset,width=ewidth,style='dashed',
                     edge_color=ecolor)
    # nx.draw_networkx_edge_labels(interest,ax=ax,pos=buses.cord,edge_labels=elabel,
    #                              font_weight='normal',font_size=25)
    
    
    ax.set_xlim(min(xpts),max(xpts))
    ax.set_ylim(min(ypts),max(ypts))
    ax.tick_params(bottom=False,left=False,labelleft=False,labelbottom=False)
    return ax



from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
axins1 = zoomed_inset_axes(ax, z, loc=2)
axins1.set_aspect(a)
axins1 = draw_subnetwork_top(axins1,G,critical,n=m)
mark_inset(ax, axins1, loc1=1, loc2=3, fc="none", ec="0.5")
axins2 = zoomed_inset_axes(ax, z, loc=3)
axins2.set_aspect(a)
axins2 = draw_subnetwork_bottom(axins2,G,critical,n=m)
mark_inset(ax, axins2, loc1=2, loc2=4, fc="none", ec="0.5")
figname = 'top-critical-edgewt-inset'+str(cont_id+1)
fig.savefig("{}{}.png".format(figpath,figname),bbox_inches='tight')

# fig=plt.figure(figsize=(50,50))
# ax=fig.add_subplot(111)
# ax = draw_subnetwork(ax,G,top_critical[0],n=3)
# figname = 'top-critical-edgewt-solo'
# fig.savefig("{}{}.png".format(figpath,figname),bbox_inches='tight')

sys.exit(0)

#%% Get edge indices
# Get edgelist and ratings
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
            
data = ''
all_edges = get_edgelist(datapath,'branchdat.csv')
for edge in top_critical:
    n_graph = get_neighbors(G,edge,n=1)
    n_edges = list(n_graph.edges(keys=True))
    crit_index = all_edges.index(edge)+1 if edge in all_edges \
        else all_edges.index((edge[1],edge[0],edge[2]))+1
    neighbor_index = [all_edges.index(e)+1 if e in all_edges \
                      else all_edges.index((e[1],e[0],e[2]))+1 for e in n_edges]
    neighbor_index.remove(crit_index)
    data += str(crit_index)+'\t'+','.join([str(x) for x in neighbor_index])+'\n'
with open(datapath+'check-criticality-'+str(K)+'.txt','w') as f:
    f.write(data)


#%% Get node indices
def get_nodelist(path,filename):
    df_nodes = pd.read_table(path+filename,sep=',')
    return df_nodes['bus_id'].tolist()

data = ''        
all_nodes = get_nodelist(datapath,'busdat.csv')
for edge in top_critical:
    n_graph = get_neighbors(G,edge,n=1)
    n_nodes = list(n_graph.nodes())
    crit_index = all_edges.index(edge)+1 if edge in all_edges \
        else all_edges.index((edge[1],edge[0],edge[2]))+1
    neighbor_index = [all_nodes.index(n)+1 for n in n_nodes]
    data += str(crit_index)+'\t'+','.join([str(x) for x in neighbor_index])+'\n'
with open(datapath+'check-load-criticality-'+str(K)+'.txt','w') as f:
    f.write(data)