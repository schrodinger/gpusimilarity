"""
These are utilities used for generating fingerprints and images.  This
implementation is RDKit-specific, using Morgan fingerprints, but it should be
trivial to change out both fingerprint and imagine generation and have all
down-stream applications use your preferred method.

NOTE:  You can speed up fingerprint creation using ipyparallel
To do this, it requires a virtualenv with ipyparallel installed
(via pip).  After sourcing the virtualenv, you can then create DBs with:
    python3 gpusim_createdb.py input.smi
YOU MUST START ipcluster from vcs-src/gpusim/python
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

# NOTE:  GPGPU backend requires fingerprint size to be sizeof(int) divisible
BITCOUNT = 1024


def add_fingerprint_bin_to_smi_line(line, trust_smiles=False):
    splitl = line.strip().split()
    try:
        smiles, cid = splitl[:2]
    except ValueError:
        raise ValueError(splitl)
    fp_binary, canon_smiles = smiles_to_fingerprint_bin(smiles,
                                            trust_smiles=trust_smiles) #noqa
    if fp_binary is None:
        return None
    return (canon_smiles, cid, fp_binary)


def split_lines_add_fp(lines, dview=None, trust_smiles=False):
    if dview is not None:
        return dview.map_sync(lambda x: add_fingerprint_bin_to_smi_line(x, trust_smiles=trust_smiles), lines)
    else:
        return map(add_fingerprint_bin_to_smi_line, lines)


def smiles_to_fingerprint_bin(smiles, trust_smiles=False):
    mol = Chem.MolFromSmiles(smiles, sanitize=(not trust_smiles))
    if mol is None:
        return None, None
    if trust_smiles:
        mol.UpdatePropertyCache()
        Chem.FastFindRings(mol)
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, BITCOUNT)

    canon_smiles = Chem.MolToSmiles(mol)
    canon_smiles = str.encode(canon_smiles)
    return DataStructs.BitVectToBinaryText(fp), canon_smiles


def smiles_to_image_file(smiles, path):
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, path)
