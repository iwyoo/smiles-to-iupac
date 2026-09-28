"""Naming benzo[cd]indole, an ortho- and peri-fused tricyclic (P-25.3.1.3):
a benzo ring shares indole's C3-C3a-C4 periphery span (letters c='3,3a',
d='3a,4', continuous per P-25.3.1.3's lettering rule), so the fused
system has one atom common to all three rings, not a simple 2-atom ortho
edge. Recognized by exact whole-molecule match, mirroring
`_azulene_parent.py`'s idiom, since only one real anchor (PubChem CID
22146199, IUPACName 'benzo[cd]indole') is available to verify a general
peri-fusion letter algorithm against.
"""

from rdkit import Chem

_BENZO_CD_INDOLE_SMILES = "C1=CC2=C3C(=C1)C=NC3=CC=C2"
_BENZO_CD_INDOLE_CANONICAL = Chem.CanonSmiles(_BENZO_CD_INDOLE_SMILES)


def has_benzo_cd_indole_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _BENZO_CD_INDOLE_CANONICAL


def name_benzo_cd_indole(mol) -> str:
    return "benzo[cd]indole"
