#!/bin/bash
# assists in the compilation of OGDF and linking
# against it in the decomposition

# build ogdf
cd OGDF
python makeMakefile.py
make

# build silp2
cd ../
make
