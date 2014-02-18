#!/usr/bin/python
'''
describes the edge graph
'''

### imports ###
import sys
import os
import logging
import networkx as nx
import numpy as np

import helpers.io as io
import optimize.orient as orient
import helpers.graphs as graphs
import helpers.misc as misc

from pygraphviz import *

### parameters ###

GRAPH_FILE = sys.argv[1]

### functions ###

def update_dists(G):
    ''' fixes distances once an orientation has been set '''
    # update distances.
    bad = 0
    to_remove = list()
    for p, q in G.edges():

        # get state and dist.
        state = G[p][q]['state']
        dist = G[p][q]['means'][state]

        # check if no dist.
        if dist == -1:
            to_remove.append((p,q))

            # re-estimate.
            width1 = G.node[p]['width']
            width2 = G.node[p]['width']
            ins_size = G[p][q]['ins_size']
            std_dev = G[p][q]['std_dev']

            # loop over distances.
            dists = np.zeros(len(G[p][q]['poses1']), dtype=np.int)
            i = 0
            for p1, p2 in zip(G[p][q]['poses1'], G[p][q]['poses2']):
                dists[i] = misc.determine_dist(p, q, state, width1, width2, p1[0], p1[1], p2[0], p2[1], ins_size)
            avg_dist = np.mean(dists)

            # remove if mean is too great.
            if avg_dist < (-1 * (ins_size + (6 * std_dev))):
                to_remove.append((p,q))
            else:
                G[p][q]['means'][state] = avg_dist

    # remove bad ones.
    G.remove_edges_from(to_remove)

    # return it.
    return G


### script ###

# load the graph.
G = nx.read_gpickle(GRAPH_FILE)

# print edge properties.
for p, q in G.edges():
    for r in G[p][q]:
        d = G[p][q][r]['dist']
        s = G[p][q][r]['state']

        print s, d




