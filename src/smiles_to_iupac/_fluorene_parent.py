"""Naming 9H-fluorene (P-25.1.2 retained name), by exact whole-molecule
match mirroring `_peri_fused_aromatic.py`'s idiom: its central ring's sp3
CH2 fails `_aromatic.py`'s all-6-membered-aromatic-ring precondition, so
no other dispatch branch ever claims it. Verified against PubChem CID 6853.
"""

from rdkit import Chem

_FLUORENE_SMILES = "C1c2ccccc2-c2ccccc21"
_FLUORENE_CANONICAL = Chem.CanonSmiles(_FLUORENE_SMILES)


def has_fluorene_parent_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _FLUORENE_CANONICAL


def name_fluorene_parent(mol) -> str:
    return "9H-fluorene"
