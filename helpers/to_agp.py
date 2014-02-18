#!/usr/bin/python
'''
translates DAG path graph to AGP.
'''

### imports ###
import sys
import os
import logging
import networkx as nx

import helpers.io as io

### definitions ###

### functions ###

def write_agp(order_file, agp_file):
	''' translates results to agp '''
	
	# load the oriented graph.
	SG = nx.read_gpickle(order_file)

	# ensure node degree is low.
	deg_list = [len(SG.neighbors(x)) for x in SG.nodes()]
	if max(deg_list) > 2:
		logging.error("is not a path")
		sys.exit(1)

	# ensure its a DAG.
	if nx.is_directed_acyclic_graph(SG) == False:
		logging.error("not a DAG?")
		sys.exit(1)

	# save it to disk.
	io.to_agp(SG, agp_file)
