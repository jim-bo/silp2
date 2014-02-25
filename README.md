silp2
=====
**Scaffolding using Integer Linear Programming**

## Overview
SILP2 is the second incarnation of the scaffolding tool developed at the University of Connecticut by James Lindsay and Dr. Ion Mandoiu, and at the Georgia State by Hamed Salooti and Dr. Alex Zelicovsky. It is built to be quick, and very flexible. The quickness comes from the usage of Nonserial Dynamic Programming which decomposes the scaffolding problem into many smaller sub-problems. The flexibility comes from the fact that it uses a generic ILP solver (CPLEX) where the constraints, objective and weights are easily tweaked.

## Usage
The tool is divided up into several sub-programs which need to be run in-order. There is a convience function to run them all. Use the "-h" argument after each of the following sub-commands for a description of their usage.
```python
python silp.py [sub-command] -h
```
1. *prep:* prepares the scaffolding input files
2. *align:* aligns the paired reads (required bowtie2) installed
3. *nodes:* creates the nodes of the scaffolding graph 
4. *edges:* adds edges to the scaffolding graph
5. *bundles:* compacts edges into bundles and computes weights
6. *decompose:* run decomposition procedure
7. *orient:* orients the contigs
8. *order:* orders the contigs
9. *gap:* computes gap sizes
10. *write_agp:* outputs the results in a common format
11. *write:* writes verbose results [debug]
12. *fasta:* writes the scaffold in fasta format N's in gaps
13. *all:* runs all the above

## Installation
This is primarily a python program, it relies on several python packages:
* numpy
* networkx
* cplex [provided by IBM](http://www-304.ibm.com/ibm/university/academic/pub/page/mem_join)

The decomposition is written in c/c++ and relies on the [OGDF](http://www.ogdf.net/doku.php) library. It must be installed and compiled. Modify the Makefile to ensure the library is included. Compilation is done by calling
```bash
make
```
in the root folder.

Bowtie2 is already required for proper alignment of reads to contigs. It must be installed and available in the $PATH variable.

## Example script.
A demonstration script file called run.sh is provided in the root directory to serve as an example on how to run the tool. A [small testcase](http://dna.engr.uconn.edu) is available to test the tool.

```bash
#!/bin/bash
# set pointer to program.
program="/path/to/silp.py"
work_dir="/path/to/data"
ref_dir="${work_dir}/ref"
asm_dir="${work_dir}/asm"
aln_dir="${work_dir}/aln"
scf_dir="${work_dir}/scf"

# align
python $program align \
       -w $scf_dir \
       -a $aln_dir \
       -p 5 \
       -c ${asm_dir}/asm.fasta \
       -s ${asm_dir}/asm.length \
       -q1 ${ref_dir}/read1.fastq \
       -q2 ${ref_dir}/read2.fastq \
       -k 2

# preprocess
python $program nodes -w ${scf_dir} -c ${asm_dir}/asm.fasta
python $program edges -w ${scf_dir} -i 3500 -s 350 -rf -s1 ${aln_dir}/read1.sam -s2 ${aln_dir}/read2.sam
python $program bundles -w ${scf_dir} -b 1 -p 90 -bup 1 -r ${aln_dir}/ant -i 3500 -s 350
python $program decompose -w ${scf_dir} -m 2500

# start time
start=$(date +%s)

# scaffold
python $program orient -w $scf_dir -z 0
python $program order -w $scf_dir
python $program gap -w $scf_dir
python $program write -w $scf_dir -a ${scf_dir}/scf.agp
python $program fasta -w $scf_dir -a ${scf_dir}/scf.agp -c ${asm_dir}/asm.fasta -f ${scf_dir}/scf.fasta

# stop time
stop=$(date +%s)
echo RUNTIME: $(expr $stop - $start) >> ${scf_dir}/scf.agp
```

## Disclaimer
This is a research tool written in a research enviroment. No support is offered and bugs may be present. Only one library size is support at this time. [CPLEX is required.](http://www-304.ibm.com/ibm/university/academic/pub/page/mem_join) a free license is available to qualified academic institutions. 

