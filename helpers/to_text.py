#!/usr/bin/python
'''
writes nodes and bundles to a text file
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

BGRAPH_FILE = sys.argv[1]
NODE_FILE = sys.argv[2]
BUNDLE_FILE = sys.argv[3]

### functions ###

### script ###

# load the bundle graph.
BG = nx.read_gpickle(BGRAPH_FILE)

# write out nodes.
with open(NODE_FILE, "wb") as fout:
	for n in BG.nodes():
		fout.write("%s\t%i\t%s\n" % (n, BG.node[n]['width'], BG.node[n]['seq']))

with open(BUNDLE_FILE, "wb") as fout:
	for p, q in BG.edges():
		fout.write("%s\t%s\t%i\t%i\t%i\t%i\n" % (p, q, BG[p][q]['bcnts'][0], BG[p][q]['bcnts'][1], BG[p][q]['bcnts'][2], BG[p][q]['bcnts'][3]))
