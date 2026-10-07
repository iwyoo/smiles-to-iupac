"""Adducts and solvates of neutral components (P-14.8, P-77.1.3(1)): the component names are joined by em dashes in
order of the seniority of classes (P-41, organic before inorganic) and the proportions follow as '(m/n)'. Components
of one class are cited alphanumerically."""

from rdkit import Chem

from ._common import UnsupportedStructure, alpha_sort_key
from ._seniority import SUFFIX_CLASS_RANK

_CLASS_PATTERNS = (
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    ("sulfonic_acid", "[SX4](=O)(=O)[OX2H1]"),
    ("ester", "[CX3](=O)[OX2][#6]"),
    ("amide", "[CX3](=O)[NX3]"),
    ("nitrile", "[CX2]#N"),
    ("aldehyde", "[CX3H1](=O)[#6]"),
    ("ketone", "[#6][CX3](=O)[#6]"),
    ("alcohol", "[#6][OX2H1]"),
    ("amine", "[NX3;!$(N[#6]=[O,S,N]);!$(N-[a])]"),
    ("ether", "[#6][OX2][#6]"),
)
_HETERO_RANK = max(SUFFIX_CLASS_RANK.values()) + 1
_HYDROCARBON_RANK = _HETERO_RANK + 1
_INORGANIC_RANK = _HYDROCARBON_RANK + 1


def _is_bare_water(frag):
    if frag.GetNumAtoms() != 1:
        return False
    atom = frag.GetAtomWithIdx(0)
    return atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 2


def _class_rank(frag):
    if not any(a.GetAtomicNum() == 6 for a in frag.GetAtoms()):
        return _INORGANIC_RANK
    ranks = [
        SUFFIX_CLASS_RANK[name]
        for name, smarts in _CLASS_PATTERNS
        if frag.HasSubstructMatch(Chem.MolFromSmarts(smarts))
    ]
    hydrocarbon = all(a.GetAtomicNum() == 6 for a in frag.GetAtoms())
    return min(ranks, default=_HYDROCARBON_RANK if hydrocarbon else _HETERO_RANK)


def _components(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) < 2 or any(a.GetFormalCharge() for a in mol.GetAtoms()):
        return None
    counts = {}
    for frag in frags:
        key = Chem.MolToSmiles(frag)
        entry = counts.setdefault(key, [frag, 0])
        entry[1] += 1
    return list(counts.items())


def has_adduct_shape(mol) -> bool:
    return _components(mol) is not None


def name_adduct(mol, namer) -> str:
    named = []
    for smiles, (frag, count) in _components(mol):
        name = "water" if _is_bare_water(frag) else namer(smiles)
        named.append((_class_rank(frag), alpha_sort_key(name), name, count))
    if len(named) == 1:
        raise UnsupportedStructure("identical components form a repeated molecule, not an adduct")
    named.sort()
    names = "—".join(name for _, _, name, _ in named)
    return f"{names} ({'/'.join(str(count) for *_, count in named)})"
