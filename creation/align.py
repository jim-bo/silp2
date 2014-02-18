#!/usr/bin/python
'''
creates node in multi scaffold graph
'''

### imports ###
import sys
import os
import logging
import subprocess
import cPickle as pickle
import networkx as nx
import numpy as np
import mmap
from operator import itemgetter
import helpers.io as io

### parameters ###

### private functions ###
def _load_lengths(path):
    """ loads lengths from file """
    
    # create dict to store info.
    lengths = dict()
    
    # loop over file.
    with open(path) as fin:
        for line in fin:
    
            # parse
            tmp = line.strip().split(" ")
            length = int(tmp[0])
            name = tmp[1]
    
            # save.
            lengths[name] = length
 
    # return info.
    return lengths
    
def _align_method(ref_idx, read_file):
    """ command line to run alignment """
    
    cmd = list()
    cmd.append("bowtie")
    cmd.append("--sam")
    cmd.append("--sam-nohead")
    cmd.append("-q")
    cmd.append("-k 20")
    cmd.append("-p 20")
    cmd.append("-k 10")
    cmd.append(ref_idx)
    cmd.append(read_file)
    return cmd

def _index_method(contig_file, index_file):
    """ command line to run alignment """
    
    cmd = list()
    cmd.append("bowtie-build")
    cmd.append(contig_file)
    cmd.append(index_file)
    return cmd

def _pair_reads(SAM1_IN_FILE, SAM2_IN_FILE, SAM1_OUT_FILE, SAM2_OUT_FILE, key_size):
    
    # memory map the SAM files1.
    fin1 = open(SAM1_IN_FILE, "r+")
    fin2 = open(SAM2_IN_FILE, "r+")

    map1 = mmap.mmap(fin1.fileno(), 0, access=mmap.ACCESS_COPY)
    map2 = mmap.mmap(fin2.fileno(), 0, access=mmap.ACCESS_COPY)

    # create lists from data.
    hitlist1 = list()
    hitlist2 = list()
    for p1, p2 in _sam_gen(map1, map2, key_size):
        hitlist1.append(p1)
        hitlist2.append(p2)

    # seek files bake to begining.
    map1.seek(0)
    map2.seek(0)

    # sort lists by name, reverse so we can pop from end.
    hitlist1.sort(key=itemgetter(1), reverse=True)
    hitlist2.sort(key=itemgetter(1), reverse=True)

    # open output files.
    fout1 = open(SAM1_OUT_FILE, "wb")
    fout2 = open(SAM2_OUT_FILE, "wb")

    # generator of pairs.
    for p1, p2 in _pair_gen(hitlist1, hitlist2):

        # load sam info from map.
        map1.seek(p1[0])
        map2.seek(p2[0])

        # write out info.
        fout1.write(map1.readline())
        fout2.write(map2.readline())

    # close output files.
    fout1.close()
    fout2.close()

    # close memmory mapped files.
    map1.close()
    map2.close()

    fin1.close()
    fin2.close()

def _sam_gen(map1, map2, KEY_SIZE):
    '''yields the SAM name and the line index'''

    # loop till end of file.
    line1 = map1.readline()
    line2 = map2.readline()
    pos1 = 0
    pos2 = 0
    while line1 != '' and line2 != '':

        # process it.
        tok1 = line1.strip().split()
        tok2 = line2.strip().split()

        # remove to key.
        if KEY_SIZE != 0:
            key1 = tok1[0][0:-KEY_SIZE]
            key2 = tok2[0][0:-KEY_SIZE]
        else:
            key1 = tok1[0]
            key2 = tok2[0]
            
        # get rname.
        rname1 = tok1[2]
        rname2 = tok2[2]

        # yield the name and line number.
        yield (pos1, key1, rname1), (pos2, key2, rname2)

        # update info.
        pos1 += len(line1)
        pos2 += len(line2)
        line1 = map1.readline()
        line2 = map2.readline()

def _pair_gen(hitlist1, hitlist2):
    ''' does an in-order walk to find pairs '''

    # loop till each list is empty.
    while len(hitlist1) > 0 and len(hitlist2) > 0:

        # peek for a match.
        if hitlist1[-1][1] == hitlist2[-1][1]:

            # check for same contig.
            if hitlist1[-1][2] != hitlist2[-1][2]:
            
                # yield it.
                yield hitlist1[-1], hitlist2[-1]

            # change left.
            hitlist1.pop()

        else:

            # pop smaller.
            if hitlist1[-1][1] < hitlist2[-1][1]:
                hitlist1.pop()
            else:
                hitlist2.pop()

### public functions ###

def create_alignment(paths, args):
    """ creates nodes
    Parameters
    ----------
    paths.node_file       : file
    args.fasta_file       : file
    """

    # local vars only.
    aln1_raw_sam = paths.tmp1_file
    aln2_raw_sam = paths.tmp2_file
    aln1_sam = args.sam1_file
    aln2_sam = args.sam2_file
    contig_file = args.contig_file
    length_file = args.length_file
    idx_dir = args.idx_dir
    idx_file = args.idx_file
    ant_file = args.ant_file
    read1_file = args.read1_file
    read2_file = args.read2_file
    ins_size = args.ins_size
    std_dev = args.std_dev
    key_size = args.key_size
    
    # create the index.
    if os.path.isdir(idx_dir) == False:
        with open('/dev/null', 'wb') as fout:
            subprocess.call(['mkdir','-p', idx_dir], stdout=fout)
            subprocess.call(_index_method(contig_file, idx_file), stdout=fout)

    # skip if necessary files are there.
    if os.path.isfile(ant_file) and os.path.isfile(aln1_sam) and os.path.isfile(aln2_sam):
        return

    # load lengths of files.
    lengths = _load_lengths(length_file)

    # create the annotation arrays.
    annotes = dict()
    for ref in lengths:
        annotes[ref] = np.zeros(lengths[ref], dtype=np.int)

    # perform alignments.
    read_paths = [read1_file, read2_file]
    out_files = [open(aln1_raw_sam, 'wb'), open(aln2_raw_sam, 'wb')]
    for i in range(2):
        
        # simplify.
        fout = out_files[i]
        
        # call alignment
        cmd = _align_method(idx_file, read_paths[i])
        
        #print ' '.join(cmd)
        #sys.exit()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        
        # setup good tracking.
        prev_qname = ""
        prev_sam = list()
        
        # iterate over results
        for line in iter(proc.stdout.readline, ''):
            
            # tokenize.
            tmp = line.strip().split("\t")
            qname = tmp[0]
            #print qname[0:-12]
            aln_status = tmp[1]
            rname = tmp[2]
            lpos = int(tmp[3])
            rpos = lpos + len(tmp[9])
            
            # skip no-align.
            if rname == "*": continue
            
            # ignore internal alignments.
            fromleft = lpos
            fromright = lengths[rname] - rpos
            if min([fromleft, fromright]) > (ins_size + (6 * std_dev)):
                continue
            
            # check if previous.
            if qname != prev_qname:
                
                # check boot condition.
                if prev_qname != "":
                    
                    # check length of previous.
                    if len(prev_sam) == 1:
                    
                        # write to file cuz its gooooood.
                        fout.write(prev_sam[0])
                        
                    else:
                        
                        # add to annotation.
                        for x in prev_sam:
                            tmp = line.strip().split()
                            qname = tmp[0]
                            rname = tmp[2]
                            lpos = int(tmp[3])
                            rpos = lpos + len(tmp[9])
                            annotes[rname][lpos:rpos] += 1
                    
                # clear prev_sam
                prev_sam = list()
                
            # update qname and same.
            prev_qname = qname
            prev_sam.append(line)

        # final check.
        if len(prev_sam) == 1:
        
            # write to file cuz its gooooood.
            fout.write(prev_sam[0])
    
        # close the file.
        fout.close()

    # clear any arrays with no entires.
    for ref in annotes:
        if np.sum(annotes[ref]) == 0:
            annotes[ref] = None

    # save the annotations.
    pickle.dump(annotes, open(ant_file, "wb" ) )
    
    # force removing from memory because it is no longer necessary.
    del annotes

    # pair the alignments.
    _pair_reads(aln1_raw_sam, aln2_raw_sam, aln1_sam, aln2_sam, key_size)

    # remove the raw file.
    subprocess.call(['rm','-f', aln1_raw_sam, aln2_raw_sam])
