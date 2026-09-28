"""Naming 4H-pyrrolo[3,2,1-ij]quinoline (julolidine's mancude parent), an
ortho- and peri-fused tricyclic (P-25.3.1.3): pyrrole shares quinoline's
N1-C8a-C8 span (letters i='8,8a', j='8a,1'), one atom (C8a) common to all
three rings. Exact whole-molecule match (PubChem CID 12311231's topology
traced against the Blue Book's `pyrrolo[3,2,1-de]acridine` bracket form).
"""

from rdkit import Chem

_PYRROLO_IJ_QUINOLINE_SMILES = "C1C=CC2=CC=CC3=C2N1C=C3"
_PYRROLO_IJ_QUINOLINE_CANONICAL = Chem.CanonSmiles(_PYRROLO_IJ_QUINOLINE_SMILES)


def has_pyrrolo_ij_quinoline_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _PYRROLO_IJ_QUINOLINE_CANONICAL


def name_pyrrolo_ij_quinoline(mol) -> str:
    return "4H-pyrrolo[3,2,1-ij]quinoline"
