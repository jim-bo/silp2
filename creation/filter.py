#!/usr/bin/python
'''
creates node in multi scaffold graph
'''

### imports ###
import sys
import os
import logging
import networkx as nx

import helpers.io as io
import helpers.misc as misc
import helpers.alignment as aln

### parameters ###

MGRAPH_FILE = sys.argv[1]

### functions ###


### script ###

# load the multi graph.
MG = nx.read_gpickle(MGRAPH_FILE)

# find 

# write to disk.
#nx.write_gpickle(MG, MGRAPH_FILE)
