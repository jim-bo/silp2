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
import optimize.order as order

### definitions ###


### functions ###

def run_ordering(paths, args):
    """ runs orientation
    Parameters
    ----------
    paths.orient_file       : file
    paths.order_file       : file
    """

    # load the bundle graph.
    BG = nx.read_gpickle(paths.bundle_file)

    # load the oriented graph.
    DG = nx.read_gpickle(paths.orient_file)

    # there are no weights in orient graph, so add them from the bundle graph
    # the weights are actually bcnts 
    for vertex1, vertex2 in DG.edges():
	DG.edge[vertex1][vertex2]["bcnts"] = BG.edge[vertex1][vertex2]["bcnts"][DG.edge[vertex1][vertex2]["state"]]

    # check it.
    for n in DG.nodes():
        if DG.node[n]['orien'] == -1:
            logging.erorr("orientatio not set")

    for p,q in DG.edges():
        if DG[p][q]['state'] == -1:
            logging.error("state not set")

    # solve the ILP.
    logging.info("solving ILP")
    ILP = order.OrderIlp("log.txt", "err.txt")
    
    # loop over each component.
    SG = nx.DiGraph()
    shorties = list()
    for subg in nx.weakly_connected_component_subgraphs(DG):
        
        # check for linear path.
        deg_list = [len(subg.neighbors(x)) for x in subg.nodes()]
        if max(deg_list) == 1 or len(subg.nodes()) < 3:
            shorties += subg.nodes()
            continue   
            
        # solve the order otherwise.
        logging.info("solving order: %s" % len(subg.nodes()))
        
        # solve it.
        ILP.load("weight", subg)
        tmp = ILP.solve()
        
        # combine it.
        SG = nx.union(SG, tmp)        
        ILP.clear()

    # solve shorties in one batch.
    subg = DG.subgraph(shorties)
    logging.info("solving shorties: %s" % len(subg.nodes()))
    ILP.load("weight", subg)
    tmp = ILP.solve()
    SG = nx.union(SG, tmp)        
    ILP.clear()
    

    # ensure node degree is low.
    logging.info("sanity check.")
    deg_list = [len(SG.neighbors(x)) for x in SG.nodes()]
    if max(deg_list) > 2:
        logging.error("is not a path")
        sys.exit(1)

    # remove cycles.
    logging.info("computing cycles")
    for subg in nx.weakly_connected_component_subgraphs(SG):
        logging.info("comp size: %i" % len(subg.nodes()))
        for cycle in nx.simple_cycles(subg):
            
            # find weakest edge.
            weakest = None
            weight = 99999
            for p, q in subg.subgraph(cycle).edges():
                s = DG[p][q]['state']
                w = BG[p][q]['bcnts'][s]
                
                if w < weight:
                    weight = w
                    weakest = p,q
                    
            # remove weakest.
            SG.remove_edge(p,q)

    # ensure its a DAG.
    if nx.is_directed_acyclic_graph(SG) == False:
        logging.error("not a DAG?")
        sys.exit(1)


    # write to disk.
    nx.write_gpickle(SG, paths.order_file)
