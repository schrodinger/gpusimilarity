# Use this Dockerfile to create DEVELOPER build.  To create a release
# build use the docker directory and follow instructions.

FROM nvidia/cuda:9.2-devel-ubuntu18.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        cuda-libraries-$CUDA_PKG_VERSION qt5-default cmake python3-pyqt5 \
        libboost-all-dev git wget python3-numpy && apt-get clean

ENV RDBASE=/rdkit-Release_2019_03_1
RUN wget https://github.com/rdkit/rdkit/archive/Release_2019_03_1.tar.gz && \
    tar zxvf Release_2019_03_1.tar.gz && cd rdkit-Release_2019_03_1 && \
    mkdir bld && cd bld && cmake -DRDK_INSTALL_INTREE=OFF \
    -DCMAKE_INSTALL_PREFIX=/usr -DPYTHON_EXECUTABLE=/usr/bin/python3 \
    -DRDK_INSTALL_STATIC_LIBS=OFF -DRDK_BUILD_CPP_TESTS=OFF .. && \
    make && make install && cd / && rm -rf Release_2019_03_1.tar.gz  \
    rdkit-Release_2019_03_1

# Non-RDKit-build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends python3-pil \
    && apt-get clean

COPY /CMake* /gpusimilarity/
COPY /*.cpp /gpusimilarity/
COPY /*.cu /gpusimilarity/
COPY /*.h /gpusimilarity/
COPY /python /gpusimilarity/python/
COPY /test /gpusimilarity/test/

RUN cd gpusimilarity && mkdir bld && cd bld && cmake ../ && make -j5
