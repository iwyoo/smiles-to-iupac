"""Linear phane names (P-26, P-52.2.5.1): four or more amplificants, two of them terminal, joined through atoms or
chains into at least seven nodes, with substituents on amplificants and bridge atoms, e.g.
2,4,6-trioxa-1,7(1),3,5(1,3)-tetrabenzenaheptaphane."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._phane_general import find_phane, is_pin_phane, name_phane_general


def has_linear_phane_shape(mol) -> bool:
    if len(Chem.GetMolFrags(mol)) != 1:
        return False
    found = find_phane(mol)
    return found is not None and not found.cyclic and is_pin_phane(found)


def name_linear_phane(mol) -> str:
    if not has_linear_phane_shape(mol):
        raise UnsupportedStructure("this structure is not a supported linear phane")
    return name_phane_general(mol)


def linear_phane_pin(mol):
    """The linear phane name of `mol`, which no semisystematic parent hydride outranks (P-52.2.5.1, P-101.1); None when
    the structure is not a supported linear phane."""
    if not has_linear_phane_shape(mol):
        return None
    try:
        return name_phane_general(mol)
    except UnsupportedStructure:
        return None
