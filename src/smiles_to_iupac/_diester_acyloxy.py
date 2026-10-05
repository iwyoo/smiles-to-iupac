"""Entry points for polyesters of one polyol (P-65.6.3.3.3): the polyol skeleton (ring system or chain) is cited as
a multivalent group before the multiplied or alphanumerically listed anions. The construction itself lives in
`_diester_ring_diyl.py` (skeleton selection, numbering) and `_diester_anions.py` (anion names). Anion names
('acetate' and the like) follow `_ester.py`.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._diester_anions import find_ester_carbons
from ._diester_ring_diyl import name_diester_ring_diyl, select_skeleton
from ._functional_prefixes import nitro_atoms

_find_ester_carbons = find_ester_carbons


def has_diester_shape(mol) -> bool:
    return len(_find_ester_carbons(mol)) >= 2


def has_polyester_of_one_polyol_shape(mol) -> bool:
    """Two or more esters with acyclic acyl carbons whose polyol has one best ring system or chain -- claimed
    ahead of the shape checks that would misread the esters."""
    nitro = nitro_atoms(mol)
    if len(Chem.GetMolFrags(mol)) > 1 or any(a.GetFormalCharge() and a.GetIdx() not in nitro for a in mol.GetAtoms()):
        return False
    matches = _find_ester_carbons(mol)
    if len(matches) < 2:
        return False
    ring_info = mol.GetRingInfo()
    if any(ring_info.NumAtomRings(acyl.GetIdx()) for acyl, _, _, _ in matches):
        return False
    try:
        return select_skeleton(mol, adjacency(mol), matches) is not None
    except UnsupportedStructure:
        return False


def name_diester_acyloxy(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    matches = _find_ester_carbons(mol)
    if len(matches) < 2:
        raise UnsupportedStructure("this module only handles two or more ester groups")
    return name_diester_ring_diyl(mol, matches)
