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
import itertools

import helpers.io as io
import optimize.orient as orient
import helpers.graphs as graphs
import helpers.misc as misc

from pygraphviz import *

### definitions ###

DECOMP_BOUND = 4

### internal functions ###

def update_dists(G):
    ''' fixes distances once an orientation has been set '''
    # update distances.
    bad = 0
    to_remove = list()
    for p, q in G.edges():

        # get state and dist.
        state = G[p][q]['state']
        dist = G[p][q]['means'][state]

        # check if no dist.
        if dist == -1:
            to_remove.append((p,q))

            # re-estimate.
            width1 = G.node[p]['width']
            width2 = G.node[p]['width']
            ins_size = G[p][q]['ins_size']
            std_dev = G[p][q]['std_dev']

            # loop over distances.
            dists = np.zeros(len(G[p][q]['poses1']), dtype=np.int)
            i = 0
            for p1, p2 in zip(G[p][q]['poses1'], G[p][q]['poses2']):
                dists[i] = misc.determine_dist(p, q, state, width1, width2, p1[0], p1[1], p2[0], p2[1], ins_size)
            avg_dist = np.mean(dists)

            # remove if mean is too great.
            if avg_dist < (-1 * (ins_size + (6 * std_dev))):
                to_remove.append((p,q))
            else:
                G[p][q]['means'][state] = avg_dist

    # remove bad ones.
    G.remove_edges_from(to_remove)

    # return it.
    return G

def _build_fixings(n, RG, SG, parents, children, comp):
    ''' build stuff to fix ILP with '''
    # build list of fixed variables.
    fixlist = list()
    if len(parents) == 1:

        # simplify.
        parent = parents[0]
        cut = RG[parent][n]['cut']

        # add fixed combos.
        if len(cut) == 1:
            fixlist.append((cut,(0,)))
        elif len(cut) > 1:
            fixlist.append((cut,(0,0)))
            fixlist.append((cut,(0,1)))

    if len(fixlist) == 0:
        fixlist.append(None)

    # build list of objective modifications.
    objdict1 = dict()
    objdict2 = dict()
    for child in children:

        # simplify.
        cut = RG[n][child]['cut']

        # handle single cut.
        if len(cut) == 1:

            # get types.
            c = list(cut).pop()
            key0 = (0,)
            key1 = (1,)
            v0 = SG.node[child]['sols'][key0]['obj']
            v1 = SG.node[child]['sols'][key1]['obj']

            # bootstrap types.
            if c not in objdict1:
                objdict1[c] = [0,0]

            # save value.
            objdict1[c][0] += v0
            objdict1[c][1] += v1

        elif len(cut) == 2:

            # get types.
            t = list(cut)
            c0 = t.pop()
            c1 = t.pop()
            key00 = (0,0)
            key01 = (0,1)
            key10 = (1,0)
            key11 = (1,1)
            v00 = SG.node[child]['sols'][key00]['obj']
            v01 = SG.node[child]['sols'][key01]['obj']
            v10 = SG.node[child]['sols'][key10]['obj']
            v11 = SG.node[child]['sols'][key11]['obj']

            # bootstrap values.
            if (c0,c1) not in objdict2:
                objdict2[(c0,c1)] = [0,0,0,0]

            # save values
            objdict2[(c0,c1)][0] += v00
            objdict2[(c0,c1)][1] += v01
            objdict2[(c0,c1)][2] += v10
            objdict2[(c0,c1)][3] += v11

    # object objective modification list.
    objlist1 = [ (c, objdict1[c][0], objdict1[c][1]) for c in objdict1]
    objlist2 = [ (c0, c1, objdict2[(c0,c1)][0], objdict2[(c0,c1)][1], objdict2[(c0,c1)][2], objdict2[(c0,c1)][3]) for c0,c1 in objdict2]

    return fixlist, objlist1, objlist2

def _blank_solution(RG, root):
    ''' creates a blank solution graph '''

    # create graph.
    SG = nx.DiGraph()

    # add root.
    SG.add_node(root, sols={():None})

    # build graph top down.
    for q in nx.dfs_preorder_nodes(RG, source=root):

        # skip root.
        if q == root: continue

        # add node.
        SG.add_node(q, sols=dict())

        # get its parent.
        p = RG.predecessors(q)[0]
        cut = RG[p][q]['cut']

        # build combos.
        fixlist = list()
        if len(cut) == 1:
            fixlist.append((0,))
            fixlist.append((1,))
        elif len(cut) > 1:
            fixlist.append((0,0))
            fixlist.append((0,1))
            fixlist.append((1,0))
            fixlist.append((1,1))

        # add each instance.
        for key in fixlist:
            SG.node[q]['sols'][key] = None

        # connect this to parents.
        SG.add_edge(p, q, cut=cut)

    # return graph.
    return SG

def _decomp_solve(BG, RG, ILP, level=0):
    ''' solves ILP using decomposition: recursive '''

    # sanity check.
    if nx.is_weakly_connected(RG) == False:
        logging.error("logic error scotty")
        sys.exit(1)

    # use topological sort to find root.
    root = nx.topological_sort(RG)[0]

    # lookup table for NSDP
    SG = _blank_solution(RG, root)

    # try to solve each node.
    for n in nx.dfs_postorder_nodes(RG, source=root):

        # simplify.
        comp = RG.node[n]['comp']
        subr = RG.node[n]['graph']
        subb = BG.subgraph(comp)
        parents = RG.predecessors(n)
        children = RG.successors(n)

        # sanity check parents.
        if len(parents) > 1:
            logging.error("decomposition isn't a tree!")
            sys.exit(1)

        # get fixings.
        fixlist, objlist1, objlist2 = _build_fixings(n, RG, SG, parents, children, comp)

        # solve each combo.
        for fix in fixlist:

            # check for direct solve or recursion.
            if len(comp) <= DECOMP_BOUND or subr == None or level == 1:

                if len(comp) > 8000:
                    logging.warning("solving a large component: %s %i" % (n, len(comp)))
                    #sys.exit()

                # solve it.
                logging.debug("solving: %s: %i: %s" % (n,len(comp), str(fix)))
                ILP.load(subb)
                ILP.fix(fix)
                
                ILP.objmod1(objlist1)
                ILP.objmod2(objlist2)
                sol = ILP.solve()
                ILP.clear()
            else:

                # recurse solve.
                logging.debug("recurse: %s: %i: %s" % (n,len(comp), str(fix)))
                sol = _decomp_solve(subb, subr, ILP, level=level+1)

            # check solution.
            for v in comp:
                if v not in sol['orien']:
                    logging.info("not all nodes solved: %s" % v)
                    sys.exit(1)

            # save solution with no cut.
            if fix == None:

                # sanity check... should be root.
                if n != root:
                    logging.error("needs to be root")
                    logging.error(str(n))
                    logging.error(str(root))
                    sys.exit(1)

                # save it to graph.
                SG.node[n]['sols'][()] = sol

            else:

                # inverse soltuion.
                invsol = {'orien':dict(), 'state':dict(), 'obj':None}
                for v in sol['orien']:
                    invsol['orien'][v] = 1 - sol['orien'][v]
                for key in sol['state']:
                    invsol['state'][key] = sol['state'][key]
                invsol['obj'] = sol['obj']

                # save key.
                key1 = fix[1]
                if len(key1) == 1:
                    key2 = (1-key1[0],)
                else:
                    key2 = (1-key1[0],1-key1[1])

                # save it to graph.
                SG.node[n]['sols'][key1] = sol
                SG.node[n]['sols'][key2] = invsol


    # top down apply solution.
    cursol = {'orien':dict(), 'state':dict(), 'obj':None}
    for p in nx.dfs_preorder_nodes(SG, source=root):

        logging.info("apply solution: %s" % p)
        #print p, RG.successors(

        # check if this is root.
        if p == root:

            # just apply
            pcomp = RG.node[p]['comp']
            sol = SG.node[p]['sols'][()]
            for v in pcomp:
                cursol['orien'][v] = sol['orien'][v]
            for key in sol['state']:
                cursol['state'][key] = sol['state'][key]
            cursol['obj'] = sol['obj']

        # loop over each child.
        for q in RG.successors(p):

            # simplify.
            qcomp = RG.node[q]['comp']
            cut = RG[p][q]['cut']

            # ensure cut is solved.
            for v in cut:
                if v not in cursol['orien']:
                    logging.error("cut node should be solved: %s" % v)
                    sys.exit()

            # find solution of cut nodes in cursol.
            if len(cut) == 1:

                # get orient of cut.
                v = list(cut).pop()
                o = cursol['orien'][v]
                sol = SG.node[q]['sols'][(o,)]

            elif len(cut) == 2:

                # get orient of cut.
                t = list(cut)
                v0 = t.pop()
                v1 = t.pop()
                o0 = cursol['orien'][v0]
                o1 = cursol['orien'][v1]
                sol = SG.node[q]['sols'][(o0,o1)]

            # save associated solution.
            for v in sol['orien']:
                cursol['orien'][v] = sol['orien'][v]
            for key in sol['state']:
                cursol['state'][key] = sol['state'][key]


            # sanity check qcomp.
            for z in qcomp:
                if z not in cursol['orien']:
                    print "HEY NOT LOADED"
                    print sol
                    print qcomp
                    sys.exit()

    # sanity check.
    totest = set(cursol['orien'].keys())
    bad = 0
    for p in RG.nodes():
        for q in RG.node[p]['comp']:
            if q not in totest:

                print "missing", p, q
                bad += 1

    if bad > 1:
        logging.error("missing nodes in solve: %s")
        sys.exit()

    # return the solution.
    return cursol

def _validate_comp_post(RG):
    ''' validates connection at a certain level '''

    # use topological sort to find root.
    root = nx.topological_sort(RG)[0]

    # try to solve each node.
    for n in nx.dfs_postorder_nodes(RG, source=root):

        # dive down.
        if RG.node[n]['graph'] != None:
            _validate_comp_post(RG.node[n]['graph'])

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



def _validate_comp_pre(RG):
    ''' validates connection at a certain level '''

    # use topological sort to find root.
    root = nx.topological_sort(RG)[0]

    # try to solve each node.
    for p in nx.dfs_preorder_nodes(RG, source=root):

        # dive down.
        if RG.node[p]['graph'] != None:
            _validate_comp_pre(RG.node[p]['graph'])

        # skip if child.
        if len(RG.successors(p)) == 0:
            continue

        # check children.
        for q in RG.successors(p):

            # get sets.
            pcomp = RG.node[p]['comp']
            qcomp = RG.node[q]['comp']

            # compute cuts.
            cutGIVEN = RG[p][q]['cut']
            cutTEST = pcomp.intersection(qcomp)

            # test cut.
            if cutGIVEN != cutTEST:
                print "bad cut"
                print cutGIVEN, cutTEST
                print n, parent
                sys.exit()

### external functions ###
def run_orientation(paths, args):
    """ runs orientation
    Parameters
    ----------
    paths.bundle_file       : file
    paths.decomp_file       : file
    """

    # load the bundle graph.
    BG = nx.read_gpickle(paths.bundle_file)
    RG = nx.read_gpickle(paths.decomp_file)

    # sanity check bundle counts.
    for p, q in BG.edges():
        x = sum(BG[p][q]['bcnts'])
        if x == 0.0:
            logging.error("zero weight bundle still present")
            sys.exit(1)

    # sanity check self edges.
    for p, q in BG.edges():
        if p == q:
            #logging.error("self edge")
            #sys.exit(1)
            BG.remove_edge(p, q)

    # annotate graph before solving.
    for p, q in BG.edges():
        BG[p][q]['state'] = -1
    for p in BG.nodes():
        BG.node[p]['orien'] = -1

    # loop over each weakly connected component.
    #ILP = orient.OrientIlp("log.txt", "err.txt", "prg.txt", "sol.txt")
    ILP = orient.OrientIlp("/dev/null", "/dev/null", "/dev/null", "/dev/null", weight_mode=args.weight_mode)
    for subcc in nx.weakly_connected_component_subgraphs(RG):

        # check its consistency.
        _validate_comp_post(subcc)
        _validate_comp_pre(subcc)

        # solve the ILP using decomposition.
        sol = _decomp_solve(BG, subcc, ILP)

        # apply the solution.
        for n in sol['orien']:
            BG.node[n]['orien'] = sol['orien'][n]
        for p,q in sol['state']:
            BG[p][q]['state'] = sol['state'][(p,q)]

    # check whole status.
    for n in BG.nodes():
        if BG.node[n]['orien'] == -1:
            logging.error("orientatio not set")

    for p,q in BG.edges():
        if BG[p][q]['state'] == -1:
            logging.error("state not set")

    # compute the directed graph.
    DG = graphs.to_directed(BG)

    # write to disk.
    nx.write_gpickle(DG, paths.orient_file)

