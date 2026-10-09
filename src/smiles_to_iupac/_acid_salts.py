"""Neutral salts of acids (P-65.6.2.1): the cations cited in alphabetical order,
then the anions, each with a multiplying prefix for repeated identical ions.
The anions are named as a whole by the substitutive engine, so partially
neutralised acids come out by the preferred method of P-65.6.2.3.1: the
anionic group is the suffix and the free acid groups are 'carboxy' prefixes."""

from rdkit import Chem

from ._common import UnsupportedStructure, alpha_sort_key
from ._multiplicative_text import enclose
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


def _is_nucleotide_anion(frag):
    from ._nucleoside_substituted import has_substituted_nucleoside_name

    return any(
        a.GetFormalCharge() == -1 and any(n.GetAtomicNum() == 15 for n in a.GetNeighbors()) for a in frag.GetAtoms()
    ) and has_substituted_nucleoside_name(frag)


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
    from ._salt import _ANION_KINDS, _POLYATOMIC_ANION_KINDS
    from .core import smiles_to_iupac

    charge = -sum(a.GetFormalCharge() for a in frag.GetAtoms())
    if not _is_nucleotide_anion(frag):
        for has_shape, namer in (*_POLYATOMIC_ANION_KINDS[:-1], *_ANION_KINDS[:-1]):
            if has_shape(frag):
                try:
                    return namer(frag), charge
                except UnsupportedStructure:
                    break
    name = smiles_to_iupac(Chem.MolToSmiles(frag))
    # the counter-ion fixes the charge of an amino acid anion (P-103.2.4.3.1)
    return (name[: -len("(1–)")] if name.endswith("(1–)") else name), charge


# 'di' + a mononuclear onium would read as the dinuclear parent (diazanium = N2H5+)
_MONONUCLEAR_ONIUM = frozenset(
    {"azanium", "phosphanium", "arsanium", "stibanium", "bismuthanium", "oxidanium", "sulfanium", "selanium", "telluranium"}
)


def _counted(entries):
    """Cited text for (name, count, simple) entries in alphabetical order with multiplying prefixes."""
    parts = []
    for name, count, simple in sorted(entries, key=lambda e: alpha_sort_key(e[0])):
        if count == 1:
            parts.append(name)
        elif simple and name.isalpha() and name not in _MONONUCLEAR_ONIUM:
            parts.append(multiplying_prefix(count) + name)
        else:
            parts.append(multiplying_prefix(count, compound=True) + enclose(name))
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
            entry = cations.setdefault(key, [found[0], found[1], 0, frag.GetNumAtoms() == 1])
        elif net < 0:
            name, charge = _anion_name(frag)
            if charge != -net:
                raise UnsupportedStructure("charge of the anion is not -net")
            entry = anions.setdefault(key, [name, charge, 0, True])
        else:
            raise UnsupportedStructure("a neutral fragment is not part of a simple salt")
        entry[2] += 1
    if not cations or not anions:
        raise UnsupportedStructure("a salt needs cations and anions")
    if sum(c * n for _, c, n, _ in cations.values()) != sum(c * n for _, c, n, _ in anions.values()):
        raise UnsupportedStructure("the charges of the ions do not balance")
    return (
        _counted([(n, k, simple) for n, _, k, simple in cations.values()])
        + " "
        + _counted([(n, k, simple) for n, _, k, simple in anions.values()])
    )
