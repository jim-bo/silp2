'''
graph functions and classes.
'''
import networkx as nx
import sys
import logging

import helpers.misc as misc
from helpers.io import load_agp


def to_directed(G):

    # make a weighted directed graph.
    DG = nx.DiGraph()

    # add the nodes.
    for n in G.nodes():
        DG.add_node(n, orien=G.node[n]['orien'], width=G.node[n]['width'])

    # add directed edges.
    for p, q in G.edges():

        # get info.
        state = G[p][q]['state']
        Sp = G.node[p]['orien']
        Sq = G.node[q]['orien']

        if state == 0:
            if Sp != Sq:
                logging.error("bad info1")
                sys.exit(1)
            if Sp == 0:
                e = p, q
            else:
                e = q, p

        elif state == 1:
            if Sp == Sq:
                logging.error("bad info2")
                sys.exit(1)
            if Sp == 0:
                e = p, q
            else:
                e = q, p

        elif state == 2:
            if Sp == Sq:
                logging.error("bad info3")
                sys.exit(1)
            if Sp == 0:
                e = q, p
            else:
                e = p, q

        elif state == 3:
            if Sp != Sq:
                logging.error("bad info4")
                sys.exit(1)
            if Sp == 0:
                e = q, p
            else:
                e = p, q

        else:
            logging.error("bad info5")
            sys.exit(1)
            
        # flip our decision based on order test.
        if p > q:
            e = e[1], e[0]
            
        # add the edge.
        DG.add_edge(e[0], e[1], state=state)

    # give back new graph.
    return DG

def agp_graph_directed(fpath):
    ''' returns agp graph '''

    # load agp array.
    agp_edges = load_agp(fpath)

    # make digraph.
    G = nx.DiGraph()

    # add nodes.
    for i in range(agp_edges.size):

        # skip contigs themselves.
        if agp_edges[i]['comp_type'] != 'W': continue

        # add node info.
        name = agp_edges[i]['comp_name']
        width = agp_edges[i]['comp_stop']
        orien = agp_edges[i]['comp_orien']

        G.add_node(agp_edges[i]['comp_name'], {'width':width, 'orien':orien})

    # add edges.
    for i in range(agp_edges.size):

        # skip contigs themselves.
        if agp_edges[i]['comp_type'] != 'N': continue

        # add sorted edges.
        ctg1 = agp_edges[i-1]['comp_name']
        ctg2 = agp_edges[i+1]['comp_name']
        o1 = G.node[ctg1]['orien']
        o2 = G.node[ctg2]['orien']
        gap = agp_edges[i]['comp_stop']
        state = misc.determine_state(ctg1, ctg2, o1, o2)

        G.add_edge(ctg1, ctg2, {'gap':gap, 'state':state})

    # done.
    return G

def agp_graph_undirected(fpath):
    ''' returns agp graph '''

    # load agp array.
    agp_edges = load_agp(fpath)

    # make digraph.
    G = nx.Graph()

    # add nodes.
    for i in range(agp_edges.size):

        # skip contigs themselves.
        if agp_edges[i]['comp_type'] != 'W': continue

        # add node info.
        name = agp_edges[i]['comp_name']
        width = agp_edges[i]['comp_stop']
        orien = agp_edges[i]['comp_orien']

        G.add_node(agp_edges[i]['comp_name'], {'width':width, 'orien':orien})

    # add edges.
    for i in range(agp_edges.size):

        # skip contigs themselves.
        if agp_edges[i]['comp_type'] != 'N': continue
        if i == agp_edges.shape[0] - 1: continue

        # add sorted edges.
        ctg1 = agp_edges[i-1]['comp_name']
        ctg2 = agp_edges[i+1]['comp_name']
        o1 = G.node[ctg1]['orien']
        o2 = G.node[ctg2]['orien']
        gap = agp_edges[i]['comp_stop']
        state = misc.determine_state(ctg1, ctg2, o1, o2)

        G.add_edge(ctg1, ctg2, {'gap':gap, 'state':state})

    # done.
    return G
