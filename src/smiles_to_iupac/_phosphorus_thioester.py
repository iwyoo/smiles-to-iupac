"""Fully esterified thio-oxoacids of phosphorus: phosphorothioates, phosphonothioates and phosphinothioates
(P-67.1.3.2, P-65.6.3.2.2) and the boron analogues borothioates, boronothioates and borinothioates (P-68.1.4.1). The
central atom bears k = 0, 1 or 2 organyl groups, for phosphorus one double-bonded chalcogen, and 3 - k single-bonded
chalcogens carrying carbon groups; at least one chalcogen is sulfur. Each ester group is cited as a
separate word with the element symbol of its chalcogen as locant ('O,O-diethyl O-phenyl phosphorothioate'), then the
anion, whose prefixes are the organyl groups and whose infix counts the sulfur atoms."""

from collections import Counter

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._phosphonic_acid import _SENIOR_ACIDS
from ._substituents import format_mononuclear_prefixes, name_branch

_PHOSPHORUS = 15
_BORON = 5
_ACID = {_PHOSPHORUS: {0: "phosphor", 1: "phosphon", 2: "phosphin"}, _BORON: {0: "boro", 1: "borono", 2: "borino"}}
_MULTIPLIER = {1: "", 2: "di"}
_INFIX = {1: "thio", 2: "dithio", 3: "trithio", 4: "tetrathio"}


def _group(mol, atom):
    """(carbon substituents, double-bonded chalcogen or None, [(single chalcogen, ester carbon)]) for a thio-oxoacid
    ester phosphorus or thio-boron-acid ester boron, else None."""
    z = atom.GetAtomicNum()
    if z not in (_PHOSPHORUS, _BORON) or atom.GetFormalCharge() or atom.IsInRing():
        return None
    if atom.GetDegree() != (4 if z == _PHOSPHORUS else 3):
        return None
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    chalcogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in (8, 16)]
    if len(carbons) + len(chalcogens) != atom.GetDegree() or len(carbons) > 2:
        return None
    double = [c for c in chalcogens if mol.GetBondBetweenAtoms(atom.GetIdx(), c.GetIdx()).GetBondTypeAsDouble() == 2.0]
    single = [c for c in chalcogens if c not in double]
    if len(double) != (1 if z == _PHOSPHORUS else 0) or (double and double[0].GetDegree() != 1) or len(single) != 3 - len(carbons):
        return None
    esters = []
    for chalcogen in single:
        rest = [n for n in chalcogen.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        if chalcogen.GetFormalCharge():
            return None
        if not rest and chalcogen.GetTotalNumHs() == 1:
            esters.append((chalcogen, None))
        elif len(rest) == 1 and rest[0].GetAtomicNum() == 6:
            esters.append((chalcogen, rest[0]))
        else:
            return None
    if all(carbon is None for _, carbon in esters):
        return None
    if not any(c.GetAtomicNum() == 16 for c in (*double, *single)):
        return None
    return carbons, double[0] if double else None, esters


def _thio_phosphorus(mol):
    return [a for a in mol.GetAtoms() if _group(mol, a) is not None]


def has_phosphorus_thioester_shape(mol) -> bool:
    return len(_thio_phosphorus(mol)) == 1


def name_phosphorus_thioester(mol) -> str:
    (phosphorus,) = _thio_phosphorus(mol)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    if sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (_PHOSPHORUS, _BORON)) != 1:
        raise UnsupportedStructure("more than one phosphorus or boron atom is not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR_ACIDS):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the phosphorus ester")
    carbons, double, esters = _group(mol, phosphorus)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    sulfurs = sum(1 for c in ((double,) if double is not None else ()) + tuple(c for c, _ in esters) if c.GetAtomicNum() == 16)

    words = Counter()
    compound = {}
    hydrogens = sum(1 for _, carbon in esters if carbon is None)
    esters = [(c, r) for c, r in esters if r is not None]
    for chalcogen, carbon in esters:
        name, is_compound = name_branch(graph, carbon.GetIdx(), chalcogen.GetIdx(), halogens, aromatic, mol=mol)
        key = (chalcogen.GetSymbol(), name)
        words[key] += 1
        compound[key] = is_compound
    cited = []
    for (symbol, name), count in sorted(words.items(), key=lambda item: (alpha_sort_key(item[0][1]), item[0][0])):
        is_compound = compound[(symbol, name)]
        group = enclose(name) if is_compound else name
        multiplier = multiplying_prefix(count, compound=is_compound) if count > 1 else ""
        cited.append(f"{','.join([symbol] * count)}-{multiplier}{group}")

    entries = [name_branch(graph, c.GetIdx(), phosphorus.GetIdx(), halogens, aromatic, mol=mol) for c in carbons]
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    stem = _ACID[phosphorus.GetAtomicNum()][len(carbons)]
    anion = f"{prefix}{stem}{'' if stem.endswith('o') else 'o'}{_INFIX[sulfurs]}ate"
    if phosphorus.GetAtomicNum() == _BORON and not hydrogens and len({c.GetSymbol() for c, _ in esters}) == 1:
        cited = [text.split("-", 1)[1] for text in cited]
    if hydrogens:
        cited.append(f"{_MULTIPLIER[hydrogens]}hydrogen")
    return " ".join([*cited, anion])
