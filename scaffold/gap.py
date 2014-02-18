#!/usr/bin/python
'''
computes distances between directed scaffold graph.
'''

### imports ###
import sys
import os
import logging
import numpy as np
import networkx as nx

import helpers.io as io
import optimize.order as order

### definitions ###


### functions ###

def compute_distance(paths, args):
    """ runs orientation
    Parameters
    ----------
    paths.order_file       : file
    paths.edge_file        : file
    paths.gap_file         : file
    """

    # load the graphs.
    BG = nx.read_gpickle(paths.bundle_file)
    EG = nx.read_gpickle(paths.edge_file)
    DG = nx.read_gpickle(paths.order_file)

    # check it.
    for n in DG.nodes():
        if DG.node[n]['orien'] == -1:
            logging.erorr("orientatio not set")

    for p,q in DG.edges():
        if DG[p][q]['state'] == -1:
            logging.error("state not set")
            
    # loop over each edge.
    missing_cnt = 0
    missing_list = list()
    for p, q in DG.edges():
        
        # simplify states.
        if DG[p][q]['state'] == 0 or DG[p][q]['state'] == 3:
            stype = 0
        else:
            stype = 1
        
        # get comparitble edges.
        glist = list()
        for e in EG[p][q]:
            
            if EG[p][q][e]['state'] == 0 or EG[p][q][e]['state'] == 3:
                ttype = 0
            else:
                ttype = 1
                
            # add edge if compatible.
            if stype == ttype:
                glist.append(e)
            
        # sanity check this.
        if len(glist) == 0:                
            glist = list()
            for e in EG[p][q]:
                glist.append(e)
            missing_cnt += 1
            missing_list.append((p,q))
            
        # compute average.
        tmp = list()
        for e in glist:
            tmp.append(EG[p][q][e]['dist'])
        avg = np.average(np.array(tmp))
            
        # save it.
        DG[p][q]['dist'] = avg

    # remove missing.
    logging.warning("removing edges with no support: %i" % missing_cnt)
    DG.remove_edges_from(missing_list)

    # write to disk.
    nx.write_gpickle(DG, paths.gap_file)    
    logging.info("missing good edges: %i" % missing_cnt)
