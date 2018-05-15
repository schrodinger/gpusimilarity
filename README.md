# gpusimilarity

## Dependencies
* RDKit (At Python level, not compilation)
* Qt 5.2+
* PyQt
* Cuda SDK, expected in default installation location /usr/local
* qmake / C++11 capable compiler

## Building
```bash
qmake && make
```

## Running
### For basic json-response http endpoint:
`python3 python/fastsim_server.py <fingerprint fsim file>`

### For testing (insecure):
`python3 python/fastsim_server.py <fingerprint fsim file> --http_interface`

### For generating databases:
Easiest from rdkit conda with pyqt installed:

```python3 python/fastsim_createdb.py <input smi.gz file> <fingerprint fsim file>```

### For debugging Cuda server, avoiding python/http server altogether:
```bash
./fastsimserver <dbname>.fsim
python3 python fastsim_search.py <dbname>
```
Note:  No .fsim extension is used for fastsim_search.py

This may be useful to determine if the backend is having Cuda/GPU problems.
