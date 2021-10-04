#!/bin/bash
bash devtools/install.sh

export PATH=`pwd`/anaconda/bin:$PATH
source activate gpu_sim

mkdir bld
cd bld
ccmake ../gpusimilarity
make -j5
ctest
