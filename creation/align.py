'''
alignment utility functions
'''
import os
import sys
import subprocess
import logging
import mmap
import gzip
import multiprocessing
from operator import itemgetter
import numpy as np


## public functions ##

def pair_alignment(paths, args):
    """ creates the alignment """

    # validate parameters.
    assert os.path.isdir(args.base_dir), 'base_dir'
    assert os.path.isfile(args.ctg_fasta), 'ctg_fasta'
    assert os.path.isfile(args.read1_sam), 'read1_fastq'
    assert os.path.isfile(args.read2_sam), 'read2_fastq'
    assert os.path.isfile(args.size_file), 'size_file'

    # key size.
    key_size = args.key_size

    # relavent files.
    base_dir = os.path.abspath(args.base_dir)
    size_file = os.path.abspath(args.size_file)
    ctg_fasta = os.path.abspath(args.ctg_fasta)
    in1_sam = os.path.abspath(args.read1_sam)
    in2_sam = os.path.abspath(args.read2_sam)

    read1_sam = os.path.abspath('%s/read1.sam' % base_dir)
    read2_sam = os.path.abspath('%s/read2.sam' % base_dir)

    names1_npy = os.path.abspath('%s/name1.npy' % base_dir)
    names2_npy = os.path.abspath('%s/name2.npy' % base_dir)
    sort1_npy = os.path.abspath('%s/sort1.npy' % base_dir)
    sort2_npy = os.path.abspath('%s/sort2.npy' % base_dir)

    ant_dir = '%s/ant' % base_dir
    idx_dir = '%s/index' % base_dir
    idx_file = '%s/index' % idx_dir

    # ensure annotation dir exists.
    subprocess.call(['mkdir', '-p', ant_dir])

    # compute name sizes.
    names_size1 = _name_size(in1_sam)
    names_size2 = _name_size(in2_sam)

    # check if sorted is present.
    if os.path.isfile(sort1_npy) == False:

        # create / load name array.
        if os.path.isfile(names1_npy) == False:
            logging.info("creating name array 1")
            names1 = _extract_names(in1_sam, names_size1, key_size)
            _names(file_name=names1_npy, data=names1)
        else:
            logging.info("loading name array 1")
            names1 = _names(file_name=names1_npy)

        # sort it.
        logging.info("sorting name array 1")
        names1.sort(order=['name'])
        _names(file_name=sort1_npy, data=names1)
        del names1
        subprocess.call(["rm", "-f", names1_npy])

    # check if sorted is present.
    if os.path.isfile(sort2_npy) == False:

        # create / load name array.
        if os.path.isfile(names2_npy) == False:
            logging.info("creating name array 2")
            names2 = _extract_names(in2_sam, names_size2, key_size)
            _names(file_name=names2_npy, data=names2)
        else:
            logging.info("loading name array 2")
            names2 = _names(file_name=names2_npy)

        # sort it.
        logging.info("sorting name array 2")
        names2.sort(order=['name'])
        _names(file_name=sort2_npy, data=names2)
        del names2
        subprocess.call(["rm", "-f", names2_npy])

    # create sizes.
    sizes = dict()
    with open(size_file, "rb") as fin:
        lines = fin.readlines()
    for line in lines:
        sz, name = line.strip().split()
        sz = int(sz)
        sizes[name] = sz

    # create the annotation arrays.
    annotes = dict()
    for ref in sizes:
        annotes[ref] = np.zeros(sizes[ref], dtype=np.int)

    # do work.
    _dual_loop(sort1_npy, sort2_npy, in1_sam, in2_sam, read1_sam, read2_sam, annotes)

    # save repeat annotation., ant_dir
    for ref in annotes:

        # create name.
        fname = '%s/%s.npy' % (ant_dir, ref)

        # look for existing.
        if os.path.isfile(fname):
            tmp = np.load(fname)
            annotes[ref] = annotes[ref] + tmp

        # serialize it.
        np.save(fname, annotes[ref])

def create_alignment(paths, args):
    """ creates the alignment """

    # validate parameters.
    assert os.path.isdir(args.base_dir), 'base_dir'
    assert os.path.isfile(args.ctg_fasta), 'ctg_fasta'
    assert os.path.isfile(args.read1_fastq), 'read1_fastq'
    assert os.path.isfile(args.read2_fastq), 'read2_fastq'

    # relavent files.
    base_dir = os.path.abspath(args.base_dir)
    ctg_fasta = os.path.abspath(args.ctg_fasta)
    read1_fastq = os.path.abspath(args.read1_fastq)
    read2_fastq = os.path.abspath(args.read2_fastq)

    tmp1_sam = os.path.abspath('%s/tmp1.sam' % base_dir)
    tmp2_sam = os.path.abspath('%s/tmp2.sam' % base_dir)

    ant_dir = '%s/ant' % base_dir
    idx_dir = '%s/index' % base_dir
    idx_file = '%s/index' % idx_dir

    # build index if not present.
    if os.path.isdir(idx_dir) == False:
        subprocess.call(["mkdir", "-p", idx_dir])
        create_idx(ctg_fasta, idx_file)

    # remove annotation dir if present.
    if os.path.isdir(ant_dir) == True:
        subprocess.call(["rm", "-rf", ant_dir])
    subprocess.call(["mkdir", "-p", ant_dir])

    # perform alignment.
    cmd1 = ['bowtie2','--reorder', '-k', '10', '-q','-p',str(args.num_cpu), '-x', idx_file, '-U', read1_fastq, '-S', tmp1_sam]
    cmd2 = ['bowtie2','--reorder', '-k', '10', '-q','-p',str(args.num_cpu), '-x', idx_file, '-U', read2_fastq, '-S', tmp2_sam]
    #print ' '.join(cmd1)
    subprocess.call(cmd1)
    #print ' '.join(cmd2)
    subprocess.call(cmd2)

def create_idx(asm_fasta, index_file):
    """     make bowtie2 index
    Parameters:
    -----------
        asm_fasta           : str
        index_file          : str
    """

    # run the command.
    subprocess.call(['bowtie2-build', '-f', asm_fasta, index_file])

## internal functions ##


def _dual_loop(sort1_npy, sort2_npy, in1_sam, in2_sam, out1_sam, out2_sam, annotes):
    """ extract unique alignments, pairs them and annotate repeats"""

    # open SAM files.
    sam1 = open(in1_sam, "rb")
    sam2 = open(in2_sam, "rb")
    out1 = open(out1_sam, "wb")
    out2 = open(out2_sam, "wb")

    # create iterators.
    itr1 = _uniq_gen(sort1_npy, sam1, annotes)
    itr2 = _uniq_gen(sort2_npy, sam2, annotes)

    # first git.
    u1 = itr1.next()
    u2 = itr2.next()
    cnt = 0
    while u1 != None and u2 != None:

        # peek for a match.
        if u1['name'] == u2['name']:

            # seek to it.
            sam1.seek(u1['row'])
            sam2.seek(u2['row'])
            out1.write(sam1.readline())
            out2.write(sam2.readline())

            # change both.
            u1 = itr1.next()
            u2 = itr2.next()

        else:

            # pop smaller.
            if u1['name'] > u2['name']:
                u2 = itr2.next()
            else:
                u1 = itr1.next()

        # die after 5
        cnt += 1
        #if cnt > 5: break

    # close them.
    sam1.close()
    sam2.close()
    out1.close()
    out2.close()


def _names(file_name=None, data=None, size=None, name_size=None):
    """ returns pointer to mapped file """

    if size != None and name_size != None:
        return np.zeros(size,  dtype=np.dtype([('name','S%d' % name_size),('row',np.int)]))
    elif file_name != None and data == None:
        return np.load(file_name)
    elif file_name != None and data != None:
        np.save(file_name, data)
    else:
        logging.error("bad usage")
        sys.exit(1)

def _name_size(file_path):
    """ guess string size """
    # determine name size.
    with open(file_path, "rb") as fin:
        for line1 in fin:
            if line1[0] == '@': continue
            name_size = len(line1.split("\t")[0]) + 10
            break

    return name_size

def _extract_names(file_name, name_size, key_size):
    """ builds numpy array of name hits"""

    # count lines.
    logging.info("reading lines")
    with open(file_name, "rb") as fin:
        size = 0
        for line in fin:
            if line[0] == '@': continue
            size += 1
            #if size > 10000000: break

    # allocate array.
    names = _names(size=size, name_size=name_size)

    # copy data into array.
    logging.info("copying data")
    with open(file_name, "rb") as fin:

        offset = 0
        idx = 0
        for line1 in fin:

            # skip header.
            if line1[0] == '@':
                offset += len(line1)
                continue

            # tokenize.
            tokens = line1.split("\t")

            # skip no map.
            if tokens[2] == "*":
                offset += len(line1)
                continue

            # operate.
            if key_size == 0:
                names[idx]['name'] = tokens[0]
            else:
                names[idx]['name'] = tokens[0][0:-key_size]
            names[idx]['row'] = offset

            # reset.
            idx += 1
            offset += len(line1)

    # resize.
    names.resize(idx)

    # return the size.
    return names



def _uniq_gen(names_npy, sam, annotes):
    """ generator for unique reads in list """

    # create mmap object.
    mmap = np.load(names_npy, mmap_mode='r')

    # setup buffered loop.
    buffstep = 10000000
    buffo = 0
    buffs = buffstep
    if buffo + buffstep > mmap.shape[0]:
        buffs = mmap.shape[0] - buffo

    # buffer loop.
    while buffo < mmap.shape[0]:

        # make buffer.
        logging.info("unique: buffering: %d %d" % (buffo, buffs))
        names = mmap[buffo:buffs]

        # iterate over non-boundry cases.
        for i in range(1, names.shape[0]-1):

            # must not match its neighbors.
            if names[i-1]['name'] != names[i]['name'] and names[i+1]['name'] != names[i]['name']:
                yield names[i]
            else:
                # annotate repeat.
                sam.seek(names[i]['row'])
                tokens = sam.readline().split("\t")
                ctg = tokens[2]
                start = int(tokens[3])
                stop = start + len(tokens[9])
                annotes[ctg][start:stop] = 1

            buffo += 1

        # check the first one.
        if names[0]['name'] != names[1]['name']:
            yield names[i]
        else:
            sam.seek(names[i]['row'])
            tokens = sam.readline().split("\t")
            ctg = tokens[2]
            start = int(tokens[3])
            stop = start + len(tokens[9])
            annotes[ctg][start:stop] = 1
        buffo += 1

        # check the last one.
        if names[-1]['name'] != names[-2]['name']:
            yield names[i]
        else:
            sam.seek(names[i]['row'])
            tokens = sam.readline().split("\t")
            ctg = tokens[2]
            start = int(tokens[3])
            stop = start + len(tokens[9])
            annotes[ctg][start:stop] = 1
        buffo += 1

        # update for buffer.
        if buffo + buffstep > mmap.shape[0]:
            buffs = buffo + (mmap.shape[0] - buffo)
        else:
            buffs = buffo + buffstep

    # yield poison pill
    yield None
