#!/usr/bin/python
'''
creates node in multi scaffold graph
'''

### imports ###
import sys
import os
import logging
import networkx as nx
import BioPyMM.sam as sam

import helpers.misc as misc

### functions ###

def create_edges(paths, args):
    """ creates edges
    Parameters
    ----------
    paths.node_file       : file
    args.sam1_file_file   : file
    args.sam2_file_file   : file
    args.pair_mode        : string
    args.ins_size         : int
    args.std_dev          : int
    paths.edge_file       : string
    """

    # load the multi graph.
    EG = nx.read_gpickle(paths.node_file)

    # add edges to the multigraph.
    for sam1, sam2, pos1, pos2 in sam.pair_gen(args.sam1_file, args.sam2_file):

        # get distance.
        p = sam1.RNAME
        q = sam2.RNAME
        width1 = EG.node[p]['width']
        width2 = EG.node[q]['width']

        # simplify.
        p = sam1.RNAME
        q = sam2.RNAME
        op, oq = misc.get_orien(sam1, sam2, args.pair_mode)
        state = misc.get_state(p, q, op, oq)
        dist = misc.get_dist(p, q, state, width1, width2, sam1.LPOS, sam1.RPOS, sam2.LPOS, sam2.RPOS, args.ins_size)

        # increment coverage.
        EG.node[p]['cov'] += sam1.RPOS - sam1.LPOS
        EG.node[q]['cov'] += sam2.RPOS - sam2.LPOS

        # skip self edges.
        if p == q:
            continue

        # add edge accordinly.
        EG.add_edge(p, q, dist=dist, state=state, left1=sam1.LPOS, right1=sam1.RPOS, left2=sam2.LPOS, right2=sam2.RPOS, ins_size=args.ins_size, std_dev=args.std_dev)

    # compute average coverage.
    for p in EG.nodes():
        EG.node[p]['cov'] = float(EG.node[p]['cov']) / float(EG.node[p]['width'])

    # ensure all edges have stuff.
    for p, q in EG.edges():
        if EG.node[p]['cov'] == 0.0 or EG.node[q]['cov'] == 0.0:
            print p, q, EG.node[p]['cov'], EG.node[q]['cov']
            for z in EG[p][q]:
                print EG[p][q][z]
            logging.error("shouldn't have happened")
            sys.exit(1)
            

    # write to disk.
    nx.write_gpickle(EG, paths.edge_file)
