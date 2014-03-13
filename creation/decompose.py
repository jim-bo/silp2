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
import subprocess
from operator import itemgetter

# hack for metis.
#os.environ['METIS_DLL'] = '/usr/local/lib/libmetis.so'
#import metis

import helpers.io as io
import optimize.orient as orient
import helpers.graphs as graphs
import helpers.misc as misc

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )

### definitions ###

DECOMP_BOUND = 500
EXT_OUT = open("/dev/null")
DECOMP_0_PROG = "%s/bin/decomp0" % '/'.join(os.path.realpath(__file__).split("/")[0:-2])
DECOMP_1_PROG = "%s/bin/decomp1" % '/'.join(os.path.realpath(__file__).split("/")[0:-2])

### functions ###

def test_bi():
    ''' creates a bi-connected component '''

    # make a biconnected component.
    G = nx.Graph()
    G.add_edge(0,1)
    G.add_edge(0,3)
    G.add_edge(0,2)
    G.add_edge(1,3)
    G.add_edge(3,2)
    G.add_edge(2,4)
    G.add_edge(4,5)
    G.add_edge(5,2)
    G.add_edge(3,7)
    G.add_edge(6,7)
    G.add_edge(6,8)
    G.add_edge(9,8)
    G.add_edge(9,7)
    return G

def test_tri():
    ''' creates a tri-connected component '''

    # make a biconnected component.
    G = nx.Graph()

    G.add_edge(1,2)
    G.add_edge(1,3)
    G.add_edge(2,4)
    G.add_edge(3,4)
    G.add_edge(5,4)
    G.add_edge(5,6)
    G.add_edge(5,7)
    G.add_edge(6,8)
    G.add_edge(7,8)
    G.add_edge(2,6)
    G.add_edge(3,7)

    G.add_edge(10,11)
    G.add_edge(10,12)
    G.add_edge(13,12)
    G.add_edge(13,11)
    G.add_edge(13,14)
    G.add_edge(15,14)
    G.add_edge(16,14)
    G.add_edge(16,17)
    G.add_edge(15,17)
    G.add_edge(11,15)
    G.add_edge(12,16)

    G.add_edge(18,19)
    G.add_edge(18,20)
    G.add_edge(21,20)
    G.add_edge(21,19)
    G.add_edge(21,22)
    G.add_edge(23,22)
    G.add_edge(24,22)
    G.add_edge(24,25)
    G.add_edge(23,25)
    G.add_edge(23,19)
    G.add_edge(20,24)

    G.add_edge(0,1)
    G.add_edge(0,10)
    G.add_edge(0,18)

    G.add_edge(9,8)
    G.add_edge(9,17)
    G.add_edge(9,25)

    return G

def make_subg(G, active):
    ''' returns subgraph '''

    # make a list of it.
    comp = list(active)

    # make new graph.
    NG = nx.Graph()

    # add nodes.
    nluf = dict()
    nlur = dict()
    for i in range(len(comp)):
        NG.add_node(i)
        nluf[comp[i]] = i
        nlur[i] = comp[i]

    # add edges.
    for p, q in G.edges():

        # skip if both not active.
        if p not in active or q not in active:
            continue

        # add the edge.
        NG.add_edge(nluf[p], nluf[q])

    # return it.
    return NG, nluf, nlur

def load_decomp(in_file, nlu, prefix):
    ''' load decomposition results '''

    # load the data.
    fin = open(in_file, "rb")
    lines = fin.readlines()
    fin.close()

    # build directed graph.
    DG = nx.DiGraph()

    # tokenize once.
    tokens = [line.strip().split() for line in lines]

    # load the nodes.
    for token in tokens:
        if token[0] != "N": continue

        # build set of ids.
        comp = frozenset([nlu[int(x)] for x in token[2::]])

        # create index of this comp.
        idx = "%s_%i" % (prefix, int(token[1]))

        # add info.
        DG.add_node(idx, comp=comp, graph=None)

    # add the directed edges.
    for token in tokens:
        if token[0] != "E": continue

        # get ids.
        s = "%s_%i" % (prefix, int(token[1]))
        t = "%s_%i" % (prefix, int(token[2]))
        cut = frozenset([nlu[int(x)] for x in token[3::]])

        # add the edge.
        DG.add_edge(s, t, cut=cut)

    # return the decomposition.
    return DG

def write_graph(G, out_file):
    ''' writes graph in simple format'''

    with open(out_file, "wb") as fout:

        # write the number of nodes and edges.
        fout.write("%i\t%i\n" % (G.number_of_nodes(), G.number_of_edges()))

        # write edges.
        for p, q in G.edges():
            fout.write("%i\t%i\n" % (p, q))

def decomp0(G, tmp1_file, tmp2_file, msize=None):
    ''' returns connected components '''

    # do decomposition.
    comps = nx.connected_components(G)

    # create the decomposition graph.
    DC = nx.DiGraph()

    # loop over connected components.
    idx = 0
    for comp in comps:

        # freeze the components.
        comp = frozenset(comp)

        # compute further decomp if necessary.
        if len(comp) > DECOMP_BOUND:
            dg = decomp1(G, comp, tmp1_file, tmp2_file, msize=msize)
        else:
            dg = None

        # add node to DC.
        DC.add_node("con_%i" % idx, comp=frozenset(comp), graph=dg)
        idx += 1

    # return the graph.
    return DC

def decomp1(G, comp, tmp1_file, tmp2_file, msize=None):
    ''' bi-connected decomposition '''

    # create active subgraph.
    subg, nluf, nlur = make_subg(G, comp)

    # serialize this to disk.
    write_graph(subg, tmp1_file)

    # execute decomposition.
    cmd = [DECOMP_0_PROG, tmp1_file, tmp2_file]
    if subprocess.call(cmd, stdout=EXT_OUT) != 0:
        logging.error("error in biconnected decomposition")
        sys.exit(1)

    # create decomposition graph.
    DC = load_decomp(tmp2_file, nlur, "bicon")

    # loop over each node in DC.
    for n in DC.nodes():

        # grab frozen component.
        comp = DC.node[n]['comp']

        # compute further decomp if necessary.
        if len(comp) > DECOMP_BOUND:

            # compute decomposition.
            dg = decomp2(G, comp, tmp1_file, tmp2_file, msize=msize)

        else:
            dg = None



        # modify node in DC.
        DC.node[n]['graph'] = dg

    # return the graph.
    return DC

def decomp2(G, comp, tmp1_file, tmp2_file, msize=None):
    ''' tri-connected decomposition '''

    # create active subgraph.
    subg, nluf, nlur = make_subg(G, comp)

    # serialize this to disk.
    write_graph(subg, tmp1_file)

    # execute decomposition.
    cmd = [DECOMP_1_PROG, tmp1_file, tmp2_file]
    if subprocess.call(cmd, stdout=EXT_OUT) != 0:
        logging.error("error in triconnected decomposition")
        sys.exit(1)

    # create decomposition graph.
    DC = load_decomp(tmp2_file, nlur, "tricon")

    # inform us of the largest component size.
    largest = -1
    lcomp = None
    for n in DC.nodes():
        if len(DC.node[n]['comp']) > largest:
            largest = len(DC.node[n]['comp'])
            lcomp = DC.node[n]['comp']
    logging.info("largest component: %d" % largest)

    # engage heuristic if necessary.
    if msize != None and len(lcomp) > msize:

        # break the graph more.
        _heuristic(G, comp)

        # note heuristic was applied.
        G.graph['redo'] = True

    # return the graph.
    return DC

def _heuristic(G, comp):
    ''' breaks component by increasing bundle size at dense cores '''

    # build subgraph.
    subg = G.subgraph(comp)
    #G.remove_edges_from(subg.edges())
    # break using balanced cuts.
    #(edgecuts, parts) = metis.part_graph(subg, 2)
    #print edgecuts
    #sys.exit()

    # rank nodes by connectivity.
    nranks = dict()
    for n in subg.nodes():
        nranks[n] = len(subg.neighbors(n))
    nranks = sorted(nranks.items(), key=itemgetter(1))

    # take top 10%
    tcut = int(len(nranks) * .05) + 1
    totrim = set([x[0] for x in nranks[-tcut::]])

    # trim any nodes with more than 5 connections.
    for n in list(totrim):
        for q in subg.neighbors(n):
            totrim.add(q)
    totrim = list(totrim)

    # build further subgraph.
    subg2 = subg.subgraph(totrim)

    # find minimum bundle size + 1
    mbs = 10000
    for p, q in subg2.edges():
        for i in range(4):
            b = subg2[p][q]['bcnts'][i]
            if b != 0 and b < mbs:
                mbs = b

    # increase minimum bundle size by 1.
    mbs += 1

    # check filter
    toremove = list()
    for p, q in subg2.edges():

        # zero out weak bundles.
        for i in range(4):
            if subg2[p][q]['bcnts'][i] < mbs:
                subg2[p][q]['bcnts'][i] = 0

        # check for edge removal.
        if sum(subg2[p][q]['bcnts']) < 1:
            toremove.append((p,q))

    # remove bad edges.
    G.remove_edges_from(toremove)
    logging.info("removed %i edges in heuristic mode" % len(toremove))

def _compact_inner(DG):
    ''' merges components if possible '''

    # recursive call to compact next level shit.
    for n in DG.nodes():
        if DG.node[n]['graph'] != None:
            _compact_inner(DG.node[n]['graph'])

    # repeat until all paths merged.
    while 1 == 1:

        # loop over each edge and check for path.
        merged = 0
        for p, q in DG.edges():

            # maybe nodes were removed already?
            if DG.has_node(p) == False or DG.has_node(q) == False:
                continue

            # don't compact node with or adj to subgraph.
            if DG.node[p]['graph'] != None or DG.node[q]['graph'] != None:
                continue

            # get childs components and kids.
            pcomp = set(DG.node[p]['comp'])
            qcomp = set(DG.node[q]['comp'])


            # size check.
            if len(pcomp) + len(qcomp) > DECOMP_BOUND:
                continue

            # merge into parent.
            DG.node[p]['comp'] = frozenset(pcomp.union(qcomp))

            # connect grandkids.
            grandkids = DG.successors(q)
            for grandkid in grandkids:
                cut = DG.node[p]['comp'].intersection(DG.node[grandkid]['comp'])
                DG.add_edge(p, grandkid, cut=cut)

            # remove node and edge come with it.
            DG.remove_node(q)
            merged += 1

        # check if we break.
        if merged == 0:
            break

def _compact_outter(DG):
    ''' merges components if possible '''

    # identify singles.
    singles = list()
    for n in DG.nodes():

        # skip connected, big and recursed.
        comp = DG.node[n]['comp']
        if len(DG.neighbors(n)) > 0: continue
        if len(comp) > DECOMP_BOUND: continue
        if DG.node[n]['graph'] != None: continue

        # note it.
        singles.append((n,len(comp)))

    # sort by small to high.
    singles = sorted(singles, key=itemgetter(1), reverse=True)

    # merge singletons.
    n, comp = singles.pop()
    curname = n
    DG.node[curname]['comp'] = set(DG.node[curname]['comp'])
    while len(singles) > 0:

        # pop next.
        n, lcomp = singles.pop()
        curset = DG.node[curname]['comp']
        nowset = set(DG.node[n]['comp'])

        # check if we can add to current.
        if len(curset) + len(nowset) < DECOMP_BOUND:

            # just concatinate.
            DG.node[curname]['comp'] = DG.node[curname]['comp'].union(nowset)

            # remove node.
            DG.remove_node(n)

        else:
            # freeze previouse.
            DG.node[curname]['comp'] = frozenset(DG.node[curname]['comp'])
            DG.node[curname]['graph'] = None

            # update current keeper.
            curname = n
            DG.node[curname]['comp'] = set(DG.node[curname]['comp'])

    # freeze last.
    DG.node[curname]['comp'] = frozenset(DG.node[curname]['comp'])
    DG.node[curname]['graph'] = None


def _validate_comp(RG):
    ''' validates connection at a certain level '''

    # use topological sort to find root.
    root = nx.topological_sort(RG)[0]

    # try to solve each node.
    for n in nx.dfs_postorder_nodes(RG, source=root):

        # dive down.
        if RG.node[n]['graph'] != None:
            _validate_comp(RG.node[n]['graph'])

        # check for parent.
        parent = RG.predecessors(n)

        # skip if root
        if len(parent) == 0:
            if n != root:
                logging.error("bad root, no cookie")
                sys.exit()
            continue
        parent = parent[0]

        # get components.
        pcomp = RG.node[parent]['comp']
        ccomp = RG.node[n]['comp']

        # compute cuts.
        cutGIVEN = RG[parent][n]['cut']
        cutTEST = pcomp.intersection(ccomp)

        # test cut.
        if cutGIVEN != cutTEST:
            print "bad cut"
            print cutGIVEN, cutTEST
            print n, parent
            sys.exit()

def decompose(paths, args):
    """ runs decomposition
    Parameters
    ----------
    paths.bundle_file       : file
    paths.tmp1_file         : file
    paths.tmp2_file         : file
    paths.decomp_file       : file
    args.msize              : integer
    """

    # load the bundle graph.
    logging.info("loading info")
    BG = nx.read_gpickle(paths.bundle_file)
    #BG = test_bi()
    #BG = test_tri()

    # run decomposition until satisfied.
    BG.graph['redo'] = False
    while 1 == 1:

        # decomposition.
        DC = decomp0(BG, paths.tmp1_file, paths.tmp2_file, msize=args.msize)

        # check if only once.
        if args.msize == None or BG.graph['redo'] == False:
            break
        elif BG.graph['redo'] == True:
            BG.graph['redo'] = False

        # remove temp files.
        if os.path.isfile(paths.tmp1_file) == True:
            subprocess.call(["rm","-f",paths.tmp1_file])
        if os.path.isfile(paths.tmp2_file) == True:
            subprocess.call(["rm","-f",paths.tmp2_file])

    # compact decomposition.
    _compact_outter(DC)
    for subcc in nx.weakly_connected_component_subgraphs(DC):

        # call recursive compaction.
        _compact_inner(DC)

    # verify decomposition.
    for subcc in nx.weakly_connected_component_subgraphs(DC):

        # check its consistency.
        _validate_comp(subcc)

    # write to disk.
    nx.write_gpickle(DC, paths.decomp_file)
    nx.write_gpickle(BG, paths.bundle_file)
