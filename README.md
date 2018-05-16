# gpusimilarity

## Dependencies
* RDKit (At Python level, not compilation)
* Qt 5.2+
* PyQt
* Cuda SDK, expected in default installation location /opt/cuda
* cmake 3.10.2+ / C++11 capable compiler
* Boost
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
`python3 ${SRC_DIR}/python/fastsim_server.py <fingerprint fsim file>`

### For testing (insecure):
From build directory:
`python3 ${SRC_DIR}/python/fastsim_server.py <fingerprint fsim file> --http_interface`

### For generating databases:
Easiest from rdkit conda with pyqt installed:

From source python directory:
```python3 fastsim_createdb.py <input smi.gz file> <fingerprint fsim file>```

### For debugging Cuda server, avoiding python/http server altogether:
```bash
From build directory:
./fastsimserver <dbname>.fsim
python3 python ${SRC_DIR}/python/fastsim_search.py <dbname>
```
Note:  No .fsim extension is used for fastsim_search.py

This may be useful to determine if the backend is having Cuda/GPU problems.
