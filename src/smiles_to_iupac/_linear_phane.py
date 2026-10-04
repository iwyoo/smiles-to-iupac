"""Linear phane names (P-26, P-52.2.5.1): four or more benzene rings, two of them terminal, joined through single
C/O/S/N nodes or chains with at least seven nodes, with substituents on rings and bridge atoms, e.g.
2,4,6-trioxa-1,7(1),3,5(1,3)-tetrabenzenaheptaphane."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._phane_general import _walks, find_phane, name_phane_general


def _shape(mol):
    found = find_phane(mol)
    if found is None:
        return None
    rings, bridges, branches, cyclic = found
    if cyclic or len(rings) < 4:
        return None
    nodes = _walks(rings, bridges, cyclic)[0]
    return found if len(nodes) >= 7 else None


def has_linear_phane_shape(mol) -> bool:
    return len(Chem.GetMolFrags(mol)) == 1 and _shape(mol) is not None


def name_linear_phane(mol) -> str:
    if _shape(mol) is None:
        raise UnsupportedStructure("this structure is not a supported linear phane")
    return name_phane_general(mol)
