"""Cyclophane names (P-26): macrocycles of benzene superatoms and bridges, delegated to the general phane engine."""

from ._phane_general import find_phane, name_phane_general

def has_cyclophane_name(mol) -> bool:
    found = find_phane(mol)
    return found is not None and found[3]


def name_cyclophane(mol) -> str:
    return name_phane_general(mol)
