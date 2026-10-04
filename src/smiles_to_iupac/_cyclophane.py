"""Cyclophane names (P-26, P-52.2.5.1): cyclic phanes with six nodes, two or more amplificants and a mancude one,
delegated to the general phane engine."""

from ._common import UnsupportedStructure
from ._phane_general import find_phane, is_fused_cyclophane, is_pin_phane, name_phane_general
from ._pin import mark


def has_cyclophane_name(mol) -> bool:
    found = find_phane(mol)
    return found is not None and found.cyclic and is_pin_phane(found)


def name_cyclophane(mol) -> str:
    return name_phane_general(mol)


def name_nonpreferred_cyclophane(mol) -> str:
    """The phane name of a cyclophane excluded from P-52.2.5.1 only by ortho or peri fusion to its ring, for structures
    with no supported fusion or bridged fused name."""
    found = find_phane(mol)
    if found is None or not is_fused_cyclophane(found):
        raise UnsupportedStructure("this structure is not a fused cyclophane")
    return mark(name_phane_general(mol), "P-52.2.5.2 prefers a fusion or bridged fused name to this phane name")
