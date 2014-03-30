#!/usr/bin/python
'''
creates bundle graph from filtered multigraph
'''

### imports ###
import sys
import os
import logging
import networkx as nx
import numpy as np
import scipy.stats as stats
import cPickle
import pickle

import helpers.io as io
import helpers.misc as misc


### definitions ###

### functions ###

def compress_edges(MG, p, q):
    ''' compresses the edges '''

    # check for types.
    bcnts = [0, 0, 0, 0]
    for z in MG[p][q]:
        bcnts[MG[p][q][z]['state']] += 1

    # build numpy arrays for each distance type.
    bdists = list()
    for i in range(4):
        bdists.append(np.zeros(bcnts[i], dtype=np.float))

    # populate array with distances.
    bidxs = [0, 0, 0, 0]
    for z in MG[p][q]:
        state = MG[p][q][z]['state']
        dist = MG[p][q][z]['dist']
        bdists[state][bidxs[state]] = dist
        bidxs[state] += 1

    # compute bundle info.
    devs = list()
    means = list()
    mins = list()
    maxs = list()
    for i in range(4):
        if bdists[i].shape[0] <= 0:
            devs.append(-1)
            means.append(-1)
            mins.append(-1)
            maxs.append(-1)
        else:
            devs.append(np.std(bdists[i]))
            means.append(np.mean(bdists[i]))
            mins.append(bdists[i].min())
            maxs.append(bdists[i].max())

    # return summaries.
    return bcnts, bdists, devs, means, mins, maxs


def _load_reps(file_path):
    ''' loads repeat info from cpickle'''

    # no weights.
    if file_path == None:
        return dict()

    # try dictionary emthod.
    if os.path.isdir(file_path) == True:
        reps = dict()
        for f in os.listdir(file_path):
            n = f.replace(".npy","")
            try:
                reps[n] = np.load("%s/%s" % (file_path, f))
            except:
                continue
        return reps

    # get weights.
    try:
        with open(file_path) as fin:
            return cPickle.load(fin)
    except:
        logging.warning("unable to load repeat pickle, ignoring weights")
        return dict()

def create_bundles(paths, args):
    """ creates bundles
    Parameters
    ----------
    paths.edge_file       : string
    args.bundle_size          : int
    args.pthresh              : int
    args.bup              : int
    """

    # load repeat annotations.
    repcnts = _load_reps(args.rep_file)

    # load the multi graph.
    MG = nx.read_gpickle(paths.edge_file)

    # create bundle graph.
    BG = nx.Graph()

    # add nodes.
    for n in MG.nodes():
        BG.add_node(n, MG.node[n])

    # build set of adjacencies.
    adjset = set()
    for p, nbrs in MG.adjacency_iter():
        for q in nbrs:
            adjset.add(tuple(sorted([p,q])))


    adjset = pickle.load(open("adjset", "rb"))

	
    # compute bundles from adjacencies.
    zerod = 0
    zcnt = 0
    ztot = len(adjset)
    for p, q in adjset:
        
        #logging.info("progress: %d of %d" % (zcnt, ztot))
        zcnt += 1

        # sanity check.
        if MG.node[p]['cov'] == 0.0 or MG.node[q]['cov'] == 0.0:
            logging.error("how can this happen?")
            sys.exit()

        # bundle size check.
        bsize = len(MG[p][q])
        if bsize < args.bundle_size:
            continue

        # group by insert size.
        groups = dict()
        std_devs = dict()
        for z in MG[p][q]:
            ins_size = MG[p][q][z]['ins_size']
            if ins_size not in groups:
                groups[ins_size] = list()
                std_devs[ins_size] = MG[p][q][z]['std_dev']
            groups[ins_size].append(z)

        # loop over groups.
        for ins_size in groups:

            # compress info.
            bcnts, bdists, devs, means, mins, maxs = compress_edges(MG, p, q)

            # compute weights.
            cov =  1 - abs(MG.node[p]['cov'] - MG.node[q]['cov']) / (MG.node[p]['cov'] + MG.node[q]['cov'])

            # swap bdists for python lists.
            for i in range(len(bdists)):
                bdists[i] = list(bdists[i])

            # add start stop info.
            poses1 = list()
            poses2 = list()
            for z in MG[p][q]:
                tmp = MG[p][q][z]
                poses1.append((tmp['left1'], tmp['right1']))
                poses2.append((tmp['left2'], tmp['right2']))

            # create bundle.
            if BG.has_edge(p, q):
                logging.error("can't have multiple insert sizes between same node")
                sys.exit(1)

            # zero out negative distances.
            avgs = [np.average(bdists[i]) for i in range(4)]
            for i in range(4):
                if avgs[i] == np.nan:
                    bcnts[i] = 0.0
                if avgs[i] < -2 * args.bundle_size:
                    bcnts[i] = 0.0
                    zerod += 1

            # don't add it if no support.
            if np.sum(bcnts) == 0:
                continue

            #BG.add_edge(p, q, bcnts=bcnts, bdists=bdists, devs=devs, means=means, mins=mins, maxs=maxs, ins_size=ins_size, std_dev=std_devs[ins_size], poses1=poses1, poses2=poses2)
            BG.add_edge(p, q, bcnts=bcnts, bdists=bdists, ins_size=ins_size, std_dev=std_devs[ins_size], cov=cov)

    # start the slimming.
    logging.info("starting repeat based slimming")

    # do repeat mods.
    track_upped = 0
    track_remed = 0
    track_ogedg = len(BG.edges())
    idxs = np.zeros(1)
    if repcnts != dict():

        # create repeat distrib.
        repavgs = np.zeros(len(repcnts), dtype=np.dtype([('name','S256'),('avg',np.float)]))
        i = 0
        for name in repcnts:

            # save the name.
            repavgs[i]['name'] = name

            # skip no repeat info.
            if name not in repcnts or repcnts[name] == None:
                repavgs[i]['avg'] = 0
                i += 1
                continue

            # take the average over ins_size + 6 (std_dev)
            d = args.ins_size + (6 * args.std_dev)
            if repcnts[name].shape[0] < d:
                repavgs[i]['avg'] = np.average(repcnts[name])
            else:
                r = range(0,d)+range(len(repcnts[name])-d,len(repcnts[name]))
                repavgs[i]['avg'] = np.average(repcnts[name][r])
            i += 1

        # compute the cutoff threshold.
        score = stats.scoreatpercentile(repavgs[:]['avg'], args.pthresh)
        idxs = repavgs[:]['avg'] > score

        # look at each bundle and see if the repeats necessitates attention.
        for p, q in BG.edges():

            # get index of pairs.
            idp = np.where(repavgs[:]['name'] == p)[0]
            idq = np.where(repavgs[:]['name'] == q)[0]

            # skip if both not high.
            if idxs[idp] == False and idxs[idq] == False:
                continue

            # get score.
            scp = repavgs[idp]['avg']
            scq = repavgs[idq]['avg']

            # check if this bundle needs attention.
            if max(scp, scq) > score:
                track_upped += 1

                # it gets its minumm bundle size upped.
                for i in range(len(BG[p][q]['bcnts'])):

                    # clear if it doesn't meet criteria.
                    if BG[p][q]['bcnts'][i] < args.bundle_size + args.bup:
                        BG[p][q]['bcnts'][i] = 0

                # remove bundle if no support.
                if np.sum(BG[p][q]['bcnts']) == 0:
                    track_remed += 1
                    BG.remove_edge(p,q)
    else:
        logging.info('no repeat information supplied')
        
    # add repeat weights.
    for p, q in BG.edges():
        
        # create weight.
        BG[p][q]['u'] = [0.0] * 4
        
        # sum weights.
        for z in MG[p][q]:
            left1 = MG[p][q][z]['left1']
            left2 = MG[p][q][z]['left2']
            right1 = MG[p][q][z]['right1']
            right2 = MG[p][q][z]['right2']            
            
            cntl = np.sum(repcnts[p][left1:left2])
            cntr = np.sum(repcnts[p][right1:right2])
            
            try:
                propl = 1.0 - (float(cntl) / float(left2-left1))
                propr = 1.0 - (float(cntr) / float(right2-right1))
            except:
                continue
                
            # add average.
            p_k = (propl + propr) / 2.0
            
            # add it.
            BG[p][q]['u'][MG[p][q][z]['state']] += p_k

    # note the modifications due to filtering.
    logging.info("contigs with repeat regions in %.2f threshold: %i of %i" % (args.pthresh, np.sum(idxs), len(idxs)))
    logging.info("bundles effected by repeats: %i of %i" % (track_upped, track_ogedg))
    logging.info("bundles removed by repeats: %i of %i" % (track_remed, track_ogedg))
    logging.info("bundles removed by neg dist: %i" % (zerod))
    logging.info("total bundles: %i" % (len(BG.edges())))

    # write to disk.
    nx.write_gpickle(BG, paths.bundle_file)
