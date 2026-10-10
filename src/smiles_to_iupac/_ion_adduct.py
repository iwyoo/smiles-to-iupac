"""Adducts that mix ions with neutral components or whose ions do not balance (P-14.8.1, P-14.8.2, P-77.1.3(1)):
balanced ions are named as the salt they form and cited as one component; ions that do not balance are cited one by one
with their charges. Components are ordered and counted as in _adduct.py."""

from functools import reduce
from math import gcd

from rdkit import Chem

from ._adduct import _MAIN_GROUP_NAMES, _class_rank, _inorganic_name, _is_bare_water
from ._common import UnsupportedStructure, alpha_sort_key

_HALIDES = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}


def _charge(frag):
    return sum(a.GetFormalCharge() for a in frag.GetAtoms())


def _components(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) < 2:
        return None
    counts = {}
    for frag in frags:
        entry = counts.setdefault(Chem.MolToSmiles(frag), [frag, 0])
        entry[1] += 1
    return [(key, frag, n) for key, (frag, n) in counts.items()]


def _ions_and_neutrals(mol):
    components = _components(mol)
    if components is None:
        return None
    ions = [c for c in components if _charge(c[1])]
    neutrals = [c for c in components if not _charge(c[1])]
    return (ions, neutrals) if ions else None


def _balanced(ions):
    return sum(_charge(frag) * n for _, frag, n in ions) == 0


def has_ion_adduct_shape(mol) -> bool:
    found = _ions_and_neutrals(mol)
    return found is not None and (bool(found[1]) or not _balanced(found[0]))


def _atom_ion_name(frag):
    from ._coordination import _METAL_NAMES

    atom = frag.GetAtomWithIdx(0)
    if frag.GetNumAtoms() != 1 or atom.GetIsotope() or atom.GetTotalNumHs():
        return None
    charge = atom.GetFormalCharge()
    if charge < 0:
        return _HALIDES.get(atom.GetAtomicNum()) if charge == -1 else None
    element = _METAL_NAMES.get(atom.GetAtomicNum()) or _MAIN_GROUP_NAMES.get(atom.GetAtomicNum())
    return f"{element}({charge}+)" if element else None


def _component_name(key, frag, namer):
    if _is_bare_water(frag):
        return "water"
    name = _inorganic_name(frag) or _atom_ion_name(frag)
    if name is not None:
        return name
    if any(a.GetNumRadicalElectrons() for a in frag.GetAtoms()):
        raise UnsupportedStructure("a radical beside ions is not supported in an adduct name")
    return namer(key)


def name_ion_adduct(mol, namer) -> str:
    found = _ions_and_neutrals(mol)
    if found is None:
        raise UnsupportedStructure("this structure holds no ion beside another component")
    ions, neutrals = found
    entries = list(neutrals)
    if _balanced(ions):
        if not neutrals:
            raise UnsupportedStructure("a balanced set of ions is a salt, not an adduct")
        formula_units = reduce(gcd, (n for _, _, n in ions))
        salt = ".".join(key for key, _, n in ions for _ in range(n // formula_units))
        salt_mol = Chem.MolFromSmiles(salt)
        if salt_mol is None:
            raise UnsupportedStructure("the ions of this structure do not form a parsable salt")
        entries.append((salt, salt_mol, formula_units))
    else:
        entries.extend(ions)
    if len(entries) < 2:
        raise UnsupportedStructure("identical components form a repeated molecule, not an adduct")
    named = []
    for key, frag, count in entries:
        name = namer(key) if len(Chem.GetMolFrags(frag)) > 1 else _component_name(key, frag, namer)
        named.append((_class_rank(frag), alpha_sort_key(name), name, count))
    named.sort(key=lambda item: item[:4])
    proportions = "/".join(str(item[3]) for item in named)
    return "—".join(item[2] for item in named) + f" ({proportions})"
