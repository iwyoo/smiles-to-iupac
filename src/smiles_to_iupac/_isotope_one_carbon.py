"""Isotopically modified one-carbon parents with a retained name (P-82.2.1, P-82.6.1.1, P-82.6.1.2): formaldehyde,
formonitrile, methanethial, methaneselenal, methanetellanal. The parent has one atom of each element, so the nuclide
descriptor carries no locant."""

from rdkit import Chem

from ._chalcogen_aldehyde import WHOLE_MOLECULES
from ._common import UnsupportedStructure
from ._isotope_labels import descriptor, split_isotopes


def _unlabelled(mol):
    split = split_isotopes(mol)
    if split is None or Chem.MolToSmiles(split[0]) not in WHOLE_MOLECULES:
        return None
    return split


def has_isotope_one_carbon_shape(mol) -> bool:
    try:
        return _unlabelled(mol) is not None
    except UnsupportedStructure:
        return False


def name_isotope_one_carbon(mol) -> str:
    clean, labels, _ = _unlabelled(mol)
    capacity = {a: clean.GetAtomWithIdx(a).GetTotalNumHs() for a in labels}
    return descriptor(labels, dict.fromkeys(labels, 1), True, capacity=capacity) + WHOLE_MOLECULES[Chem.MolToSmiles(clean)]
