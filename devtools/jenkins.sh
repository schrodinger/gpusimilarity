#!/bin/bash
bash devtools/install.sh

export PATH=`pwd`/anaconda/bin:$PATH
source activate gpu_sim

mkdir bld
cd bld
ccmake ../
make -j5
ctest
