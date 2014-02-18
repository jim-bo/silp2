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

### parameters ###

### functions ###

def create_nodes(paths, args):
    """ creates nodes
    Parameters
    ----------
    paths.node_file       : file
    args.fasta_file       : file
    """

    # read in fasta to dictionary.
    seqs = io.load_fasta(args.contig_file)

    # create graph.
    G = nx.MultiGraph()

    # add nodes to graph.
    for name, seq in seqs.items():

        # skip split names.
        tmp = name.split(" ")
        name = tmp[0]

        # add node.
        G.add_node(name, {'seq':seq, 'width':len(seq), 'cov':0})

    # write to disk.
    nx.write_gpickle(G, paths.node_file)
