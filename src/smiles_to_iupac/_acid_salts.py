"""Neutral salts of acids (P-65.6.2.1): the cations cited in alphabetical order,
then the anions, each with a multiplying prefix for repeated identical ions.
The anions are named as a whole by the substitutive engine, so partially
neutralised acids come out by the preferred method of P-65.6.2.3.1: the
anionic group is the suffix and the free acid groups are 'carboxy' prefixes."""

from rdkit import Chem

from ._acid_groups import acid_group_at
from ._carbonic_family import carbonic_anion_center
from ._common import UnsupportedStructure, alpha_sort_key
from ._numerals import multiplying_prefix

_EXTRA_CATIONS = {
    ("Ge", 4): "germanium",
    ("Sb", 3): "antimony",
    ("Bi", 3): "bismuth",
    ("Zn", 2): "zinc",
    ("Cd", 2): "cadmium",
    ("Ga", 3): "gallium",
    ("In", 3): "indium",
}


def _has_acid_anion(frag):
    for atom in frag.GetAtoms():
        if atom.GetFormalCharge() == -1 and atom.GetAtomicNum() in (8, 16, 34, 52):
            for n in atom.GetNeighbors():
                if n.GetAtomicNum() in (6, 16, 34, 52) and not n.GetFormalCharge():
                    group = acid_group_at(frag, n.GetIdx())
                    if group is not None and group.spec.anion or carbonic_anion_center(frag, n.GetIdx()):
                        return True
    return False


def _cation(frag):
    from ._salt import _cation as legacy_cation

    found = legacy_cation(frag)
    if found is not None:
        return found
    if frag.GetNumAtoms() == 1:
        atom = frag.GetAtomWithIdx(0)
        name = _EXTRA_CATIONS.get((atom.GetSymbol(), atom.GetFormalCharge()))
        if name is not None:
            return name, atom.GetFormalCharge()
    return None


def _anion_name(frag):
    from ._salt import _polyatomic_anion
    from .core import smiles_to_iupac

    polyatomic = _polyatomic_anion(frag)
    if polyatomic is not None:
        namer, magnitude = polyatomic
        return namer(frag), magnitude
    charge = -sum(a.GetFormalCharge() for a in frag.GetAtoms())
    return smiles_to_iupac(Chem.MolToSmiles(frag)), charge


def _simple(name):
    return name.isalpha()


def _counted(entries):
    """Cited text for (name, count) entries in alphabetical order with multiplying prefixes."""
    parts = []
    for name, count in sorted(entries, key=lambda e: alpha_sort_key(e[0])):
        if count == 1:
            parts.append(name)
        elif _simple(name):
            parts.append(multiplying_prefix(count) + name)
        else:
            parts.append(multiplying_prefix(count, compound=True) + f"({name})")
    return " ".join(parts)


def name_acid_salt(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) < 2:
        raise UnsupportedStructure("not a salt")
    cations, anions = {}, {}
    for frag in frags:
        net = sum(a.GetFormalCharge() for a in frag.GetAtoms())
        key = Chem.MolToSmiles(frag)
        if net > 0:
            found = _cation(frag)
            if found is None or found[1] != net:
                raise UnsupportedStructure("this cation is not named by the salt engine")
            entry = cations.setdefault(key, [found[0], found[1], 0])
        elif net < 0:
            if not _has_acid_anion(frag) and _anion_name_is_unsafe(frag):
                raise UnsupportedStructure("this anion is not an acid anion")
            name, charge = _anion_name(frag)
            if charge != -net:
                raise UnsupportedStructure("charge of the anion is not -net")
            entry = anions.setdefault(key, [name, charge, 0])
        else:
            raise UnsupportedStructure("a neutral fragment is not part of a simple salt")
        entry[2] += 1
    if not cations or not anions:
        raise UnsupportedStructure("a salt needs cations and anions")
    if not any(_has_acid_anion(Chem.MolFromSmiles(k)) for k in anions):
        raise UnsupportedStructure("no acid anion in this salt")
    if sum(c * n for _, c, n in cations.values()) != sum(c * n for _, c, n in anions.values()):
        raise UnsupportedStructure("the charges of the ions do not balance")
    return _counted([(n, k) for n, _, k in cations.values()]) + " " + _counted([(n, k) for n, _, k in anions.values()])


def _anion_name_is_unsafe(frag):
    from ._salt import _polyatomic_anion

    return _polyatomic_anion(frag) is None
