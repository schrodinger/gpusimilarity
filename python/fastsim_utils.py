"""
These are utilities used for generating fingerprints and images.  This
implementation is RDKit-specific, using Morgan fingerprints, but it should be
trivial to change out both fingerprint and imagine generation and have all
down-stream applications use your preferred method.

NOTE:  You can speed up fingerprint creation using ipyparallel
To do this, it requires a virtualenv with ipyparallel installed
(via pip).  After sourcing the virtualenv, you can then create DBs with:
    python3 fastsim_createdb.py input.smi
YOU MUST START ipcluster from vcs-src/fastsim/python
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

# NOTE:  GPGPU backend requires fingerprint size to be sizeof(int) divisible
BITCOUNT = 1024
RETURN_COUNT = 10  # must match server return count, defined in types.h


def add_fingerprint_bin_to_smi_line(line, trustSmiles=False):
    splitl = line.strip().split()
    try:
        smiles, cid = splitl[:2]
    except ValueError:
        raise ValueError(splitl)
    fp_binary = smiles_to_fingerprint_bin(smiles, trustSmiles=trustSmiles)
    if fp_binary is None:
        return None
    return (smiles, cid, fp_binary)


def split_lines_add_fp(lines, dview=None, trustSmiles=False):
    if dview is not None:
        return dview.map_sync(lambda x: add_fingerprint_bin_to_smi_line(x, trustSmiles=trustSmiles), lines)
    else:
        return map(add_fingerprint_bin_to_smi_line, lines)


def smiles_to_fingerprint_bin(smiles, trustSmiles=False):
    mol = Chem.MolFromSmiles(smiles, sanitize=(not trustSmiles))
    if mol is None:
        return None
    if trustSmiles:
        mol.UpdatePropertyCache()
        Chem.FastFindRings(mol)
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, BITCOUNT)

    return DataStructs.BitVectToBinaryText(fp)


def smiles_to_image_file(smiles, path):
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, path)
