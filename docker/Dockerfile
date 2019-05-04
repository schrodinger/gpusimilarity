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

# Use specific hash/version in checkout to ensure cacheing works
RUN git clone https://github.com/schrodinger/gpusimilarity.git && \
    cd gpusimilarity && \
    git checkout 1570f5c46ff21a6ae8238a28ca0fa3a13d8e87ac && \
    mkdir bld && cd bld && cmake ../ && make -j5

