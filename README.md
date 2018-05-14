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
`./fastsimserver <fingerprint fsim file>`

### For testing (insecure):
`./fastsimserver <fingerprint_file> --http_interface`

### For generating databases:
Easiest from rdkit conda with pyqt installed:

```python3 python/fastsim_createdb.py <input smi.gz file> <fingerprint fsim file>```
