FROM nvidia/cuda:9.2-devel-ubuntu18.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        cuda-libraries-$CUDA_PKG_VERSION qt5-default cmake python3-pyqt5 \
        libboost-all-dev git wget python3-numpy && apt-get clean

ENV RDBASE=/rdkit-Release_2019_03_1
RUN wget https://github.com/rdkit/rdkit/archive/Release_2019_03_1.tar.gz && \
    tar zxvf Release_2019_03_1.tar.gz && cd rdkit-Release_2019_03_1 && \
    mkdir bld && cd bld && cmake -DRDK_INSTALL_INTREE=OFF \
    -DCMAKE_INSTALL_PREFIX=/usr -DPYTHON_EXECUTABLE=/usr/bin/python3 .. && \
    make && make install && cd / && rm -rf Release_2019_03_1.tar.gz  \
    rdkit-Release_2019_03_1

# IMPORTANT:  Change this to remote branch you want to build
ENV BRANCH=docker

# Used to prevent github cache from stopping rebuild when branch changes
ADD https://api.github.com/repos/schrodinger/gpusimilarity/git/refs/heads/$BRANCH version.txt

RUN git clone https://github.com/schrodinger/gpusimilarity.git && \
    cd gpusimilarity && \
    git checkout $BRANCH && \
    mkdir bld && cd bld && cmake ../ && make -j5

