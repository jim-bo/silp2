#!/usr/bin/python
'''
creates node in multi scaffold graph
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

### definitions ###

DECOMP_BOUND = 50

### functions ###

def flip_existing(G):
    ''' flips all already oriented nodes '''
    for n in G.nodes():
        if G.node[n]['orien'] != -1:
            G.node[n]['orien'] = 1 - G.node[n]['orien']

def debug_orientation(paths, args):
    ''' runs the orientation '''

    # simplify.
    bundle_file = paths.bundle_file

    # load the bundle graph.
    BG = nx.read_gpickle(bundle_file)

    # sanity check self edges.
    for p, q in BG.edges():
        if p == q:
            BG.remove_edge(p, q)

    # annotate graph before solving.
    for n in BG.nodes():
        BG.node[n]['orien'] = -1
    for p, q in BG.edges():
        BG[p][q]['state'] = -1

    # remove branching nodes.
    for n in BG.nodes():
        if len(BG.neighbors(n)) > 2:
            for q in BG.neighbors(n):
                BG.remove_edge(n,q)

    # make edge assignment by connected componnent.
    in_it = False
    for comp in nx.connected_components(BG):

        # handle singelton.
        if len(comp) == 1:
            BG.node[comp[0]]['orien'] = 0
            continue

        # make subgraph.
        subg = BG.subgraph(comp)

        # find hanging node.
        for n in subg.nodes():
            if subg.neighbors(n) == 1:
                break

        if len(subg.edges()) == 2:
            for p, q in subg.edges():
                print p, q, BG[p][q]['bcnts'], p<q
            in_it = True

        # assign orientation of first.
        BG.node[n]['orien'] = 0

        # generate ordered edges.
        for p, q in nx.dfs_edges(subg, n):

            # query maximum state.
            state = np.argmax(BG[p][q]['bcnts'])
            BG[p][q]['state'] = state

            # assign orientation.
            if state == 0:
                if BG.node[p]['orien'] == 0:
                    BG.node[q]['orien'] = 0
                else:
                    BG.node[q]['orien'] = 1

            elif state == 1:
                if BG.node[p]['orien'] == 0:
                    BG.node[q]['orien'] = 1
                else:
                    BG.node[q]['orien'] = 0

            elif state == 2:
                if BG.node[p]['orien'] == 1:
                    BG.node[q]['orien'] = 0
                else:
                    BG.node[q]['orien'] = 1

            elif state == 3:
                if BG.node[p]['orien'] == 1:
                    BG.node[q]['orien'] = 1
                else:
                    BG.node[q]['orien'] = 0

        if in_it == True:
            for n in subg.nodes():
                print n, BG.node[n]['orien']
            for p, q in subg.edges():
                print p, q, BG[p][q]['state']

            print "--"
            DG = graphs.to_directed(BG.subgraph(comp))
            for p, q in DG.edges():
                print p, q
            sys.exit()

    sys.exit()

    # compute the directed graph.
    #DG = graphs.to_directed(BG)

    # write to disk.
    #nx.write_gpickle(DG, paths.orient_file)

