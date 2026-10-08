"""Condensed guanidines H2N-[C(=NH)-NH]n-H (P-66.4.1.2.2): the diamides of imidodicarbonimidic acid (n = 2),
diimidotricarbonimidic acid (n = 3) and triimidotetracarbonimidic acid (n = 4), and for n = 5 and higher the
skeletal replacement ('a') names 'triimino-tetraazanonane-diimidamide'."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, group_substituents
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
    return order, bridges


def has_condensed_guanidine_shape(mol) -> bool:
    return _chain(mol) is not None


def _terminal_nitrogens(mol, carbon, bridges):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(carbon).GetNeighbors() if n.GetIdx() not in bridges]


def _unit_entries(mol, graph, carbon, bridges, halogens, aromatic_atoms):
    """[(locant letter, name, compound)] of the substituents on the terminal nitrogens of one end unit."""
    entries = []
    for nitrogen in _terminal_nitrogens(mol, carbon, bridges):
        imino = mol.GetBondBetweenAtoms(carbon, nitrogen).GetBondTypeAsDouble() == 2.0
        for n in graph[nitrogen]:
            if n == carbon:
                continue
            name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol)
            entries.append(("N'" if imino else "N", name, compound))
    return entries


def name_condensed_guanidine(mol) -> str:
    found = _chain(mol)
    if found is None or len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("not a condensed guanidine")
    order, bridges = found
    count = len(order)
    if mol.GetNumAtoms() == 3 * count + 1:
        return _plain_name(count)
    if count != 2:
        raise UnsupportedStructure("substituted condensed guanidines beyond the biguanide are not supported yet")
    return _substituted_biguanide(mol, order, bridges)


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


def _substituted_biguanide(mol, order, bridges):
    from ._common import halogen_substituents

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    units = []
    for carbon in order:
        terminals = _terminal_nitrogens(mol, carbon, bridges)
        if len(terminals) != 2:
            raise UnsupportedStructure("an unusual condensed guanidine tautomer is not supported yet")
        if any(
            mol.GetBondBetweenAtoms(carbon, b).GetBondTypeAsDouble() == 2.0 for b in bridges if b in graph[carbon]
        ):
            raise UnsupportedStructure("the double bond to the bridging nitrogen leaves the locants ambiguous")
        units.append(_unit_entries(mol, graph, carbon, bridges, halogens, aromatic_atoms))
    ranked = sorted(units, key=lambda entries: (-len(entries), sorted(name for _, name, _ in entries)))
    positions = {}
    for digit, entries in enumerate(ranked, start=1):
        for letter, name, compound in entries:
            positions.setdefault(name, {"locants": [], "compound": compound})["locants"].append(f"{letter}{digit}")
    prefix = format_substituent_prefixes(positions) if positions else ""
    return f"{prefix}{_ACID_STEMS[2]} diamide"
