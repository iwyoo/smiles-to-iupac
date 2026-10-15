"""Condensed guanidines H2N-[C(=NH)-NH]n-H (P-66.4.1.2.2): the diamides of imidodicarbonimidic acid (n = 2),
diimidotricarbonimidic acid (n = 3) and triimidotetracarbonimidic acid (n = 4), and for n = 5 and higher the
skeletal replacement ('a') names 'triimino-tetraazanonane-diimidamide'."""

from rdkit import Chem

from ._chalcogenourea import has_group_senior_to_urea, unprimed_first
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._numerals import alkane_name, multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_ACID_STEMS = {2: "imidodicarbonimidic", 3: "diimidotricarbonimidic", 4: "triimidotetracarbonimidic"}


def _guanidine_carbons(mol):
    return [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6
        and a.GetDegree() == 3
        and not a.GetFormalCharge()
        and all(n.GetAtomicNum() == 7 for n in a.GetNeighbors())
        and sorted(mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in a.GetNeighbors()) == [1.0, 1.0, 2.0]
    ]


def _chain(mol):
    """The guanidine carbons in order along the NH bridges, or None."""
    carbons = set(_guanidine_carbons(mol))
    if len(carbons) < 2:
        return None
    bridges = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge():
            continue
        ends = [n.GetIdx() for n in atom.GetNeighbors()]
        if all(e in carbons for e in ends):
            bridges[atom.GetIdx()] = ends
    if len(bridges) != len(carbons) - 1:
        return None
    degree = {c: 0 for c in carbons}
    for ends in bridges.values():
        for e in ends:
            degree[e] += 1
    if any(d > 2 for d in degree.values()) or sorted(degree.values()).count(1) != 2:
        return None
    order = [next(c for c, d in degree.items() if d == 1)]
    used = set()
    while len(order) < len(carbons):
        step = next(
            ((b, ends) for b, ends in bridges.items() if b not in used and order[-1] in ends), None
        )
        if step is None:
            return None
        used.add(step[0])
        order.append(next(e for e in step[1] if e != order[-1]))
    chain_atoms = set(order) | set(bridges) | {n for c in order for n in _terminal_nitrogens(mol, c, bridges)}
    if any(mol.GetAtomWithIdx(a).IsInRing() for a in chain_atoms):
        return None
    return order, bridges


def has_condensed_guanidine_shape(mol) -> bool:
    return _chain(mol) is not None


def _terminal_nitrogens(mol, carbon, bridges):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(carbon).GetNeighbors() if n.GetIdx() not in bridges]


def _unit_entries(mol, graph, carbon, terminals, halogens, aromatic_atoms):
    """[(prime, name, compound)] of the substituents on the terminal nitrogens of one guanidine carbon. The =NH group of an
    end unit is N' and its -NH2 group N; where the double bond runs to a bridging nitrogen the substituted group takes the
    unprimed locant, as in the tautomers of guanidine (P-66.4.1.2.1.2)."""
    named = {
        nitrogen: [
            name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol) for n in graph[nitrogen] if n != carbon
        ]
        for nitrogen in terminals
    }
    if len(terminals) == 1:
        return [("", name, compound) for name, compound in named[terminals[0]]]
    imino = [n for n in terminals if mol.GetBondBetweenAtoms(carbon, n).GetBondTypeAsDouble() == 2.0]
    if imino:
        primed = imino[0]
    else:
        first, second = terminals
        primed = second if unprimed_first(named[first], named[second]) else first
    return [("'" if n == primed else "", name, compound) for n in terminals for name, compound in named[n]]


def name_condensed_guanidine(mol) -> str:
    found = _chain(mol)
    if found is None or len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("not a condensed guanidine")
    order, bridges = found
    count = len(order)
    if mol.GetNumAtoms() == 3 * count + 1:
        return _plain_name(count)
    return _substituted_name(mol, order, bridges)


def _plain_name(count):
    if count in _ACID_STEMS:
        return f"{_ACID_STEMS[count]} diamide"
    length = 2 * count - 1
    imino = list(range(3, length - 1, 2))
    aza = list(range(2, length, 2))
    return (
        f"{','.join(map(str, imino))}-{multiplying_prefix(len(imino))}imino-"
        f"{','.join(map(str, aza))}-{multiplying_prefix(len(aza))}aza{alkane_name(length)}-1,{length}-diimidamide"
    )


def _substituted_name(mol, order, bridges):
    """The unit at position 2k-1 of the chain carries N and N' of that position (P-66.4.1.2.2); the end units have two
    nitrogens to substitute, a middle unit one. The numbering direction gives the substituents the lowest locants, then
    the alphanumerically first prefix the lowest locant."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    core = set(order) | set(bridges) | {n for carbon in order for n in _terminal_nitrogens(mol, carbon, bridges)}
    if has_group_senior_to_urea(mol, core):
        raise UnsupportedStructure("a group senior to the condensed guanidine is the parent, not an N-substituent")
    units = []
    for index, carbon in enumerate(order):
        terminals = _terminal_nitrogens(mol, carbon, bridges)
        if len(terminals) != (2 if index in (0, len(order) - 1) else 1):
            raise UnsupportedStructure("an unusual condensed guanidine tautomer is not supported yet")
        units.append(_unit_entries(mol, graph, carbon, terminals, halogens, aromatic_atoms))
    best = None
    for direction in (units, units[::-1]):
        located = sorted(
            ((2 * index + 1, prime), name, compound)
            for index, entries in enumerate(direction)
            for prime, name, compound in entries
        )
        key = ([place for place, _, _ in located], [alpha_sort_key(name) for _, name, _ in located])
        if best is None or key < best[0]:
            best = (key, located)
    positions = {}
    for (digit, prime), name, compound in best[1]:
        positions.setdefault(name, {"locants": [], "compound": compound})["locants"].append(f"N{prime}{digit}")
    prefix = format_substituent_prefixes(positions) if positions else ""
    return f"{prefix}{_ACID_STEMS[len(order)]} diamide"
