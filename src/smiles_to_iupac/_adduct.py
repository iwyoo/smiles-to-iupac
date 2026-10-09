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
_DONOR_ATOMS = {7, 8, 15, 16, 33, 34, 51, 52}
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


def _donor_text(base, donor):
    """(donor symbol or N locant, component name citing its locants) for a donor atom that shares its element with other
    atoms of the component: the nitrogens of a urea or guanidine take the locants N, N', N'' (P-66.1.1.4.1.1,
    P-66.4.1.2.1.1)."""
    from ._guanidine import nitrogen_locants as guanidine_locants
    from ._urea import nitrogen_locants as urea_locants

    if donor.GetAtomicNum() != 7:
        return None
    for locants in (urea_locants, guanidine_locants):
        found = locants(base)
        if found is not None and donor.GetIdx() in found[0]:
            return found[0][donor.GetIdx()], found[1]
    return None


def _attachment(base, acid):
    """(donor text, acceptor symbol, component name or None) for a drawn donor-acceptor pair, cited when the donor
    component has several possible donor atoms (P-68.1.6.2); None otherwise. The donor text is the element symbol when
    that is unambiguous in the component, else the locant of the donor atom."""
    donors = [a for a in base.GetAtoms() if a.HasProp("_adduct_donor")]
    acceptors = [a for a in acid.GetAtoms() if a.HasProp("_adduct_acceptor")]
    if len(donors) != 1 or len(acceptors) != 1:
        return None
    donor, acceptor = donors[0], acceptors[0]
    candidates = [a.GetAtomicNum() for a in base.GetAtoms() if a.GetAtomicNum() in _DONOR_ATOMS]
    if len(candidates) < 2:
        return None
    if [a.GetAtomicNum() for a in acid.GetAtoms()].count(acceptor.GetAtomicNum()) != 1:
        return None
    if candidates.count(donor.GetAtomicNum()) == 1:
        return donor.GetSymbol(), acceptor.GetSymbol(), None
    located = _donor_text(base, donor)
    return None if located is None else (located[0], acceptor.GetSymbol(), located[1])


def has_adduct_shape(mol) -> bool:
    return _components(mol) is not None


def name_adduct(mol, namer) -> str:
    named = []
    for smiles, (frag, count) in _components(mol):
        name = "water" if _is_bare_water(frag) else namer(smiles)
        named.append((_class_rank(frag), alpha_sort_key(name), name, count, frag))
    if len(named) == 1:
        raise UnsupportedStructure("identical components form a repeated molecule, not an adduct")
    named.sort(key=lambda item: item[:4])
    proportions = f"({'/'.join(str(item[3]) for item in named)})"
    if len(named) == 2 and all(item[3] == 1 for item in named):
        pair = _attachment(named[0][4], named[1][4])
        if pair is not None:
            return f"{pair[2] or named[0][2]}({pair[0]}\u2014{pair[1]}){named[1][2]} {proportions}"
    joined = "\u2014".join(item[2] for item in named)
    return f"{joined} {proportions}"
