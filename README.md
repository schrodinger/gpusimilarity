# gpusimilarity

A brute-force GPU implementation of chemical fingerprint similarity searching.  Its intended use is to be kept alive as a service with an entire library loaded into graphics card memory.  It has python scripts included which use RDKit to generate fingerprints, but the C++/Cuda backend are agnostic to the data once it's been created.

## Initial Benchmark

Basic benchmarks for searching a 17M [Zinc-based](http://zinc.docking.org/) library:

Full similarity comparison against all 17M compounds with 1024bit fingerprints, including sort of results:

Tesla V100:  0.025 seconds (~680M a second)

GeForce 1080Ti:  0.05 seconds (~350M a second)

## Example integration

Here is a video of this backend being utilized for immediate-response searching inside Schrödinger's LiveDesign application:

[![GPUSimilarity Gadget](http://img.youtube.com/vi/DZhknAXXEo4/0.jpg)](http://www.youtube.com/watch?v=DZhknAXXEo4)

## Dependencies
* RDKit (At Python level, not compilation)
* Qt 5.2+ (including QtNetwork)
* PyQt
* Cuda SDK, CUDACXX env variable pointing to nvcc
* cmake 3.10.2+
* C++11 capable compiler
* Boost test libraries
* Optional: Doxygen for generating documents

## Building with CMake and running unit tests with CTest
```
From parent directory of source:
mkdir bld
cd bld
ccmake ../gpusimilarity
make -j5
ctest
```
If Cuda, boost or doxygen are not found, start ccmake with the following
options:
```
ccmake -DCMAKE_CUDA_COMPILER=/path/to/nvcc -DBOOST_ROOT=/path/to/boost/directory -DDOXYGEN_EXECUTABLE=/path/to/doxygen
```
### Generate the documentation
Install doxygen on system
```
make doc_doxygen
```
The result is in bld/doc/html

## Running
### For basic json-response http endpoint:
From build directory:
`python3 ${SRC_DIR}/python/gpusim_server.py <fingerprint fsim file>`

### For testing (insecure):
From build directory:
`python3 ${SRC_DIR}/python/gpusim_server.py <fingerprint fsim file> --http_interface`

### For generating databases:
Easiest from rdkit conda with pyqt installed:

From source python directory:
```python3 gpusim_createdb.py <input smi.gz file> <fingerprint fsim file>```

### For debugging Cuda server, avoiding python/http server altogether:
```bash
From build directory:
./gpusimserver <dbname>.fsim
python3 python ${SRC_DIR}/python/gpusim_search.py <dbname>
```
Note:  No .fsim extension is used for gpusim_search.py

This may be useful to determine if the backend is having Cuda/GPU problems.
