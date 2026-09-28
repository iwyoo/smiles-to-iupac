"""Naming azulene (P-25.1.1 item 17, a retained name for the 5-7 fused
bicyclic hydrocarbon), by exact whole-molecule match mirroring
`_fluorene_parent.py`'s idiom: its 5- and 7-membered rings both fail
`_aromatic.py`'s all-6-membered-ring precondition. Verified against
PubChem CID 9231; fusion-atom locants 3a/8a confirmed at P-25.3.3.1.2(c).
"""

from rdkit import Chem

_AZULENE_SMILES = "C1=CC2=CC=CC=CC2=C1"
_AZULENE_CANONICAL = Chem.CanonSmiles(_AZULENE_SMILES)


def has_azulene_parent_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _AZULENE_CANONICAL


def name_azulene_parent(mol) -> str:
    return "azulene"
