#!/usr/bin/python
'''
test bundle graph against an AGP for valid adjacencies.
'''

### imports ###
import sys
import os
import logging
import networkx as nx
import numpy as np

import helpers.io as io
import helpers.misc as misc
import helpers.graphs as gphs


### parameters ###

BGRAPH_FILE = sys.argv[1]
AGP_FILE = sys.argv[2]
INSERT_SIZE = int(sys.argv[3])
STD_DEV = int(sys.argv[4])


### functions ###

def find_adj(G):
	'''uses DFS from every node to find all nodes within reach.'''
	
	# make a new graph to track edges.
	UG = nx.Graph()
	
	# loop over each node.
	for n in G.nodes():
		
		# DFS from this node.
		stack = [n]
		visited = set()
		
		while len(stack) > 0:
			
			# peek at stack.
			p = stack[-1]
			
			# add edge to this.
			if UG.has_edge(n, p) == False:
				UG.add_edge(n, p)
			
			# mark as visited.
			visited.add(p)
			
			# loop over neighbors.
			hit = False
			for q in G.neighbors(p):
				
				# check if we visited this already.
				if q not in visited:
					
					
					# find shortest path.
					path = nx.shortest_path(G, n, q)
					size = 0
					for z in path[1:-1]:
						size += G.node[z]['width']
					for i in range(0, len(path[1:-1])-1):
						size += G[path[i]][path[i+1]]['gap']
						
					# see if we can reach it via dist.
					if size <= INSERT_SIZE + (3 * STD_DEV):
						stack.append(q)
						hit = True
			
			# continue if a hit.
			if hit == True:
				continue
				
			# pop it from stack.
			stack.pop()
				
	# return graph.
	return UG

### script ###

# load the graphs.
BG = nx.read_gpickle(BGRAPH_FILE)
AG = gphs.agp_graph_undirected(AGP_FILE)

# remove large edges from AGP.
to_remove = list()
for p, q in AG.edges():
	if AG[p][q]['gap'] > INSERT_SIZE + (3 * STD_DEV):
		to_remove.append((p,q))
AG.remove_edges_from(to_remove)

# build transitive adjacencies.
UG = find_adj(AG)

# turn edges into sets.
BGset = set([tuple(sorted([p,q])) for p,q in BG.edges()])
AGset = set([tuple(sorted([p,q])) for p,q in AG.edges()])
UGset = set([tuple(sorted([p,q])) for p,q in UG.edges()])

# count adjacencies.
basic = len(AGset.intersection(BGset))
trans = len(UGset.intersection(BGset))

# report.
print "Actual Edges: %i" % len(AG.edges())
print "Chosen Edges: %i" % len(BG.edges())
print "Basic Matching: %i" % basic
print "Transitive Matching: %i" % trans
