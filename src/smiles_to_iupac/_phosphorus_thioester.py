"""Fully esterified thio-oxoacids of phosphorus: phosphorothioates, phosphonothioates and phosphinothioates
(P-67.1.3.2, P-65.6.3.2.2) and the boron analogues borothioates, boronothioates and borinothioates (P-68.1.4.1). The
central atom bears k = 0, 1 or 2 organyl groups, for phosphorus one double-bonded chalcogen, and 3 - k single-bonded
chalcogens carrying carbon groups; at least one chalcogen is sulfur. Each ester group is cited as a
separate word with the element symbol of its chalcogen as locant ('O,O-diethyl O-phenyl phosphorothioate'), then the
anion, whose prefixes are the organyl groups and whose infix counts the sulfur atoms."""

from collections import Counter

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._acid_lexicon import _peroxo_word
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._phosphonic_acid import _SENIOR_ACIDS
from ._substituents import format_mononuclear_prefixes, name_branch

_PHOSPHORUS = 15
_BORON = 5
_ACID = {_PHOSPHORUS: {0: "phosphor", 1: "phosphon", 2: "phosphin"}, _BORON: {0: "boro", 1: "borono", 2: "borino"}}
_MULTIPLIER = {1: "", 2: "di"}
_INFIX = {1: "thio", 2: "dithio", 3: "trithio", 4: "tetrathio"}


def _peroxy_link(mol, link, bridge):
    """The carbon of an -O-O-R, -O-S-R, -S-O-R or -S-S-R group whose first chalcogen is `bridge`, else None."""
    if link.GetAtomicNum() not in (8, 16) or link.GetDegree() != 2 or link.GetFormalCharge():
        return None
    carbons = [n for n in link.GetNeighbors() if n.GetIdx() != bridge.GetIdx() and n.GetAtomicNum() == 6]
    return carbons[0] if len(carbons) == 1 else None


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
            esters.append((chalcogen, None, None))
        elif len(rest) == 1 and rest[0].GetAtomicNum() == 6:
            esters.append((chalcogen, rest[0], None))
        elif len(rest) == 1 and _peroxy_link(mol, rest[0], chalcogen):
            esters.append((chalcogen, _peroxy_link(mol, rest[0], chalcogen), rest[0]))
        else:
            return None
    if all(carbon is None for _, carbon, _ in esters):
        return None
    if not any(link for _, _, link in esters) and not any(c.GetAtomicNum() == 16 for c in (*double, *single)):
        return None
    return carbons, double[0] if double else None, esters


def _thio_phosphorus(mol):
    return [a for a in mol.GetAtoms() if _group(mol, a) is not None]


def _linked_pair(mol):
    """The two phosphonic thioester phosphorus atoms of a multiplicative ester, else None: each carries one carbon (the
    linker) and no free acid hydroxy group."""
    groups = _thio_phosphorus(mol)
    if len(groups) != 2 or sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (_PHOSPHORUS, _BORON)) != 2:
        return None
    if any(a.GetAtomicNum() != _PHOSPHORUS for a in groups):
        return None
    found = [_group(mol, a) for a in groups]
    if any(len(carbons) != 1 or any(carbon is None or link is not None for _, carbon, link in esters) for carbons, _, esters in found):
        return None
    return groups, found


def has_phosphorus_thioester_shape(mol) -> bool:
    return len(_thio_phosphorus(mol)) == 1 or _linked_pair(mol) is not None


def _name_linked_pair(mol, groups, found):
    from ._radical_poly import _benzene, _component

    if len(Chem.GetMolFrags(mol)) > 1 or any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("a multi-fragment or isotopically modified structure is not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR_ACIDS):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the phosphorus ester")
    sulfur_counts = {
        sum(1 for c in (*(() if double is None else (double,)), *(c for c, _, _ in esters)) if c.GetAtomicNum() == 16)
        for _, double, esters in found
    }
    if len(sulfur_counts) != 1 or 0 in sulfur_counts:
        raise UnsupportedStructure("the two acid groups differ in their chalcogen replacement")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    attachments = [(carbons[0].GetIdx(), atom.GetIdx(), 1) for atom, (carbons, _, _) in zip(groups, found)]
    part = _component(mol, attachments)
    if part is None:
        raise UnsupportedStructure("the linker of the two phosphorus acid groups is not named")
    best = None
    for order in ((0, 1), (1, 0)):
        words = {}
        for rank, index in enumerate(order):
            for chalcogen, carbon, _ in found[index][2]:
                name, is_compound = name_branch(graph, carbon.GetIdx(), chalcogen.GetIdx(), halogens, aromatic, mol=mol)
                words.setdefault(name, [is_compound, []])[1].append((rank, chalcogen.GetSymbol()))
        cited = []
        for name, (is_compound, locants) in sorted(words.items(), key=lambda item: alpha_sort_key(item[0])):
            text = ",".join(symbol + "\u2032" * rank for rank, symbol in sorted(locants))
            group = enclose(name) if is_compound else name
            multiplier = multiplying_prefix(len(locants), compound=is_compound) if len(locants) > 1 else ""
            cited.append(f"{text}-{multiplier}{group}")
        key = tuple(cited)
        if best is None or key < best:
            best = key
    infix = _INFIX[next(iter(sulfur_counts))]
    anion = f"phosphono{infix}ate"
    text = _benzene(part.text)
    linker = enclose(text) if part.has_prefix or part.has_locants else text
    return " ".join([*best, f"P,P\u2032-{linker}bis({anion})"])


def name_phosphorus_thioester(mol) -> str:
    pair = _linked_pair(mol)
    if pair is not None:
        return _name_linked_pair(mol, *pair)
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
    pairs = [(c, link) for c, _, link in esters if link is not None]
    if len(pairs) > 1:
        raise UnsupportedStructure("several peroxy ester groups on one acid centre are not supported yet")
    sulfurs = sum(
        1
        for c in ((double,) if double is not None else ()) + tuple(c for c, _, link in esters if link is None)
        if c.GetAtomicNum() == 16
    )

    words = Counter()
    compound = {}
    hydrogens = sum(1 for _, carbon, _ in esters if carbon is None)
    esters = [entry for entry in esters if entry[1] is not None]
    for chalcogen, carbon, link in esters:
        name, is_compound = name_branch(
            graph, carbon.GetIdx(), (link or chalcogen).GetIdx(), halogens, aromatic, mol=mol
        )
        key = (chalcogen.GetSymbol() + (link.GetSymbol() if link else ""), name)
        words[key] += 1
        compound[key] = is_compound
    unlocated = phosphorus.GetAtomicNum() == _BORON and not hydrogens and len(set("".join(key[0] for key in words))) == 1
    if unlocated:
        merged = Counter()
        for (_, name), count in words.items():
            merged[("", name)] += count
        words = merged
    cited = []
    for (symbol, name), count in sorted(words.items(), key=lambda item: (alpha_sort_key(item[0][1]), item[0][0])):
        is_compound = compound.get((symbol, name), any(v for (_, n), v in compound.items() if n == name))
        group = enclose(name) if is_compound else name
        multiplier = multiplying_prefix(count, compound=is_compound) if count > 1 else ""
        text = f"{multiplier}{group}"
        cited.append(text if unlocated else f"{','.join([symbol] * count)}-{text}")

    entries = [name_branch(graph, c.GetIdx(), phosphorus.GetIdx(), halogens, aromatic, mol=mol) for c in carbons]
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    stem = _ACID[phosphorus.GetAtomicNum()][len(carbons)]
    infix = _INFIX[sulfurs] if sulfurs else ""
    if pairs:
        peroxo = _peroxo_word((pairs[0][0].GetSymbol(), pairs[0][1].GetSymbol()))
        infix += f"({peroxo}ate)" if peroxo != "peroxo" else "peroxoate"
        anion = f"{prefix}{stem}{'' if stem.endswith('o') else 'o'}{infix}"
        anion += "" if infix.endswith(("ate)", "peroxoate")) else "ate"
    else:
        anion = f"{prefix}{stem}{'' if stem.endswith('o') else 'o'}{infix}ate"
    if hydrogens:
        cited.append(f"{_MULTIPLIER[hydrogens]}hydrogen")
    return " ".join([*cited, anion])


def name_boron_peroxy_ester(mol) -> str:
    """A boron acid ester with a peroxy ester group: 'O-ethyl OS-methyl phenylborono(thioperoxoate)' (P-68.1.4.1)."""
    groups = _thio_phosphorus(mol)
    if len(groups) != 1 or groups[0].GetAtomicNum() != _BORON or not any(link for _, _, link in _group(mol, groups[0])[2]):
        raise UnsupportedStructure("a boron acid ester with a peroxy ester group is required")
    return name_phosphorus_thioester(mol)
