#!/usr/bin/python
'''
The main script for scaffold algorithm SILP.
'''
### imports ###

# system
import subprocess
import warnings
import argparse
import logging
import time
import sys
import os

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s', )

# local
import creation.align
import creation.nodes
import creation.edges
import creation.bundles
import creation.decompose
import scaffold.part1
import scaffold.part2
import scaffold.gap
import helpers.io

# hack to silence argparser.
warnings.filterwarnings('ignore', category=DeprecationWarning)

### classes ###

class WorkPaths(object):
    ''' file paths for this scaffolding '''


    def __init__(self, args):
        ''' creates and validates paths'''

        # save them.
        # check working dir.
        if os.path.isdir(args.work_dir) == False:
            subprocess.call(['mkdir','-p',args.work_dir])
        try:
            self.work_dir = args.work_dir
        except AttributeError:
            self.work_dir = os.path.abspath(self.work_dir)
            
        try:
            self.contig_file = os.path.abspath(args.contig_files)
        except AttributeError:
            self.contig_file = None
            
        try:
            self.sam1_file = os.path.abspath(args.sam1_files)
            self.sam2_file = os.path.abspath(args.sam2_files)
        except AttributeError:
            self.sam1_file = None
            self.sam2_file = None

        # validate work dir.
        if self.work_dir == None or os.path.isdir(self.work_dir) == False:
            logging.error('bad work dir %s' % self.work_dir)
            sys.exit(1)
        self.work_dir = os.path.abspath(self.work_dir)

        # validate supplied and aobfiles.
        for key, val in vars(self).items():
            if key != 'work_dir':
                if val != None:
                    if os.path.isfile(val) == False:
                        logging.error('bad supplied file: %s' % val)
                        sys.exit(1)
                    else:
                        self.__dict__[key] = os.path.abspath(val)

        # user has no control over these.
        self.node_file = '%s/node.cpickle' % self.work_dir
        self.edge_file = '%s/edge.cpickle' % self.work_dir
        self.bundle_file = '%s/bundle.cpickle' % self.work_dir
        self.decomp_file = '%s/decomp.cpickle' % self.work_dir
        self.orient_file = '%s/orient.cpickle' % self.work_dir
        self.order_file = '%s/order.cpickle' % self.work_dir
        self.gap_file = '%s/gap.cpickle' % self.work_dir
        self.tmp1_file = '%s/tmp1.txt' % self.work_dir
        self.tmp2_file = '%s/tmp2.txt' % self.work_dir
        

### functions ###

def allatonce(paths, args):
    ''' runs everything with timing '''
    # prepare the folder.
    prepare_experiment(paths, args)
    
    # create the alignments.
    creation.align.create_alignment(paths, args)
    
    # create the scaffolding graph.
    creation.nodes.create_nodes(paths, args)
    creation.edges.create_edges(paths, args)
    creation.bundles.create_bundles(paths, args)
    creation.decompose.decompose(paths, args)

    tstart = time.time()
    scaffold.part1.run_orientation(paths, args)
    scaffold.part2.run_ordering(paths, args)
    scaffold.gap.compute_distance(paths, args)
    tstop = time.time()

    #helpers.io.write_agp(paths, args, runtime=tstop-tstart)
    helpers.io.write_agp(paths, args)

def prepare_experiment(paths, args):
    ''' symbolic links everything '''
        
    # symbolic link files.
    files = [args.contig_file, args.length_file, args.read1_file, args.read2_file]
    names = ["contigs.fa","contigs.length","read1.sam","read2.sam"]
    for x, y in zip(files, names):
        x = os.path.abspath(x)
        y = '%s/%s' % (os.path.abspath(args.work_dir), y)
        if os.path.isfile(y) == False:
            subprocess.call(['ln','-s',x,y])

### script ###

if __name__ == '__main__':

    # mode parser.
    main_p = argparse.ArgumentParser()
    subp = main_p.add_subparsers(help='sub-command help')

    # prepare the scaffolding.
    node_p = subp.add_parser('prep', help='prepares scaffolding experiment')
    node_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    node_p.add_argument('-c', dest='contig_file', required=True, help='contig file')
    node_p.add_argument('-l', dest='length_file', required=True, help='length file')
    node_p.add_argument('-r1', dest='read1_file', required=True, help='read1 file')
    node_p.add_argument('-r2', dest='read2_file', required=True, help='read2 file')
    node_p.set_defaults(func=prepare_experiment)    

    # create the alignment.
    aln_p = subp.add_parser('align', help='aligns reads to contigs and sets up necessary files')
    aln_p.add_argument('-w', dest='work_dir', required=True, help='scaffolding directory')
    aln_p.add_argument('-a', dest='base_dir', required=True, help='alignment directory')
    aln_p.add_argument('-p', dest='num_cpu', type=int, required=True, help='number of threads for bowtie2')
    aln_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    aln_p.add_argument('-q1', dest='read1_fastq', required=True, help='read first file')
    aln_p.add_argument('-q2', dest='read2_fastq', required=True, help='read second file')
    aln_p.set_defaults(func=creation.align.create_alignment)

    # pair an existing alignment.
    aln_p = subp.add_parser('pair', help='aligns reads to contigs and sets up necessary files')
    aln_p.add_argument('-w', dest='work_dir', required=True, help='scaffolding directory')
    aln_p.add_argument('-a', dest='base_dir', required=True, help='alignment directory')
    aln_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    aln_p.add_argument('-s1', dest='read1_sam', required=True, help='read first file')
    aln_p.add_argument('-s2', dest='read2_sam', required=True, help='read second file')
    aln_p.add_argument('-l', dest='size_file', required=True, help='size file')
    aln_p.add_argument('-k', dest='key_size', type=int, required=True, help='size of PE key at end of each read')
    aln_p.set_defaults(func=creation.align.pair_alignment)

    # node sub-parser.
    node_p = subp.add_parser('nodes', help='creates node graph from contig file')
    node_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    node_p.add_argument('-c', dest='contig_file', required=True, help='contig file')
    node_p.set_defaults(func=creation.nodes.create_nodes)

    # edge sub-parser.
    edge_p = subp.add_parser('edges', help='creates edge graph from nodes and SAM files')
    edge_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    edge_p.add_argument('-i', dest='ins_size', type=int, action='store', required=True, help='insert size')
    edge_p.add_argument('-s', dest='std_dev', type=int, action='store', required=True, help='standard deviation')
    edge_g = edge_p.add_mutually_exclusive_group(required=True)
    edge_g.add_argument('-ff', dest='pair_mode', action='store_const', const=0, help='SOLiD style -> ->')
    edge_g.add_argument('-fr', dest='pair_mode', action='store_const', const=1, help='innie style -> <-')
    edge_g.add_argument('-rf', dest='pair_mode', action='store_const', const=2, help='outtie style <- ->')
    edge_p.add_argument('-s1', dest='sam1_file', required=True, help='first file in SAM pair')
    edge_p.add_argument('-s2', dest='sam2_file', required=True, help='second file in SAM pair')
    edge_p.set_defaults(func=creation.edges.create_edges)
    
    # bundle sub-parser.
    bundle_p = subp.add_parser('bundles', help='creates bundles from edge graph')
    bundle_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    bundle_p.add_argument('-b', dest='bundle_size', type=int, action='store', required=True, help='bundle size')
    bundle_p.add_argument('-i', dest='ins_size', type=int, action='store', required=True, help='insert size')
    bundle_p.add_argument('-s', dest='std_dev', type=int, action='store', required=True, help='standard deviation')
    bundle_p.add_argument('-r', dest='rep_file', action='store', default=None, help='repeat count dir')
    bundle_p.add_argument('-p', dest='pthresh', default=90, type=int, action='store', help='percentile threshold')
    bundle_p.add_argument('-bup', dest='bup', default=1, type=int, action='store', help='up suspicious bundles by this')
    bundle_p.set_defaults(func=creation.bundles.create_bundles)

    # decomposition sub-parser.
    decomp_p = subp.add_parser('decompose', help='decomposes the graph into smaller parts')
    decomp_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    decomp_p.add_argument('-m', dest='msize', default=1, type=int, action='store', help='engage heuristic on large components.')
    decomp_p.set_defaults(func=creation.decompose.decompose)

    # orient sub-parser.
    orient_p = subp.add_parser('orient', help='orient the bundle graph')
    orient_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    orient_p.add_argument('-z', dest='weight_mode', type=int, default=0, help='weight mode')
    orient_p.set_defaults(func=scaffold.part1.run_orientation)

    # order sub-parser.
    order_p = subp.add_parser('order', help='order the directed scaffold graph')
    order_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    order_p.set_defaults(func=scaffold.part2.run_ordering)
    
    # calculate gaps.
    order_p = subp.add_parser('gap', help='estimates gaps')
    order_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    order_p.set_defaults(func=scaffold.gap.compute_distance)
    
    # write gap information for hamed
    order_p = subp.add_parser('write_gap', help='write info for computing gaps')
    order_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    order_p.set_defaults(func=helpers.io.write_gap_info)

    # write sub-parser.
    write_p = subp.add_parser('write', help='write the AGP file')
    write_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    write_p.add_argument('-a', dest='agp_file', required=True, help='agp file')
    write_p.set_defaults(func=helpers.io.write_agp)

    # write sub-parser.
    write_p = subp.add_parser('fasta', help='write the fasta file')
    write_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    write_p.add_argument('-c', dest='ctg_file', required=True, help='fastafile')
    write_p.add_argument('-a', dest='agp_file', required=True, help='agp file')
    write_p.add_argument('-f', dest='scf_file', required=True, help='fastafile')
    write_p.set_defaults(func=helpers.io.write_fasta)


    # all sub-parser.
    all_p = subp.add_parser('all', help='runs entire scaffolding operation')
    all_p.add_argument('-w', dest='work_dir', required=True, help='working directory')
    all_p.add_argument('-c', dest='contig_file', required=True, help='contig file')
    all_p.add_argument('-l', dest='length_file', required=True, help='length file')
    all_p.add_argument('-r1', dest='read1_file', required=True, help='read1 file')
    all_p.add_argument('-r2', dest='read2_file', required=True, help='read2 file')
    all_p.add_argument('-i', dest='ins_size', type=int, action='store', required=True, help='insert size')
    all_p.add_argument('-s', dest='std_dev', type=int, action='store', required=True, help='standard deviation')
    all_p.add_argument('-key', dest='key_size', type=int, required=True, help='how much to chop off end of read name')
    all_p.add_argument('-p', dest='pthresh', default=90, type=int, action='store', help='percentile threshold')
    all_p.add_argument('-bup', dest='bup', default=1, type=int, action='store', help='up suspicious bundles by this')
    all_p.add_argument('-b', dest='bundle_size', type=int, action='store', required=True, help='bundle size')
    all_g = all_p.add_mutually_exclusive_group(required=True)
    all_g.add_argument('-ff', dest='pair_mode', action='store_const', const=0, help='SOLiD style -> ->')
    all_g.add_argument('-fr', dest='pair_mode', action='store_const', const=1, help='innie style -> <-')
    all_g.add_argument('-rf', dest='pair_mode', action='store_const', const=2, help='outtie style <- ->')
    all_p.add_argument('-a', dest='agp_file', required=True, help='agp file')
    all_p.set_defaults(func=allatonce)

    # parse args.
    args = main_p.parse_args()
    args.func(WorkPaths(args), args)
