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
* cplex [provided by IBM]

However the decomposition is written in c/c++ and relies on the [OGDF](http://www.ogdf.net/doku.php) library. It must be installed and compiled. Modify the Makefile to ensure the libraries is included. Compilation is done by calling
```bash
make
```
in the root folder.



## Disclaimer
This is a research tool written in a research enviroment. No support is offered and bugs may be present. Only one library size is support at this time. CPLEX is required. 

