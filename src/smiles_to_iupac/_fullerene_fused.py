"""Organic rings and ring systems ortho fused to a fullerene or modified fullerene, cited as fusion prefixes with primed
locants in front of the fullerene parent (P-27.6.1). The fusion bond stays single and the fusion atoms take no exo double
bond, so the ring component is judged without that bond: its remaining saturated positions take indicated hydrogen
when odd in number and hydro prefixes otherwise (P-31.1.4.2.4)."""

from dataclasses import dataclass

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure
from ._fused_numbering import _locant_key
from ._fullerene_modified import modified_candidate
from ._fusion_components import _RETAINED_PREFIX, identify, skeleton
from ._numerals import multiplying_prefix

_PRIMES = ("", "′", "″", "‴", "⁗")
_DONORS = frozenset({"O", "S", "Se", "Te"})
_MAX_PIECE = 40
_MAX_ATOMS = 160


@dataclass
class _Ring:
    prefix: str
    bracketed: str
    cite: bool
    compound: bool
    options: list


def name_fullerene_derivative(mol):
    """Name of `mol` as a modified fullerene, with the fused rings it carries, else None. A cage that needs 'cyclo' or
    'seco' only to take in the atoms of a ring is named with that ring as fusion prefix instead."""
    plain = modified_candidate(mol)
    if plain is not None and not (plain.edit.seco or plain.edit.cyclo):
        return plain.text()
    fused = _fused_name(mol)
    if fused is not None:
        return fused
    return None if plain is None else plain.text()


def _fused_name(mol):
    ring_info = mol.GetRingInfo()
    if mol.GetNumAtoms() > _MAX_ATOMS or ring_info.NumRings() < 11 or sum(len(r) == 5 for r in ring_info.AtomRings()) < 5:
        return None
    graph = nx.Graph((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds())
    pieces = _pieces(mol, graph)
    if not pieces:
        return None
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return None
    rings = []
    for a, b, atoms in pieces:
        ring = _ring(mol, kekule, graph, a, b, atoms)
        if ring is None:
            return None
        rings.append(ring)
    fusion_atoms = {x for a, b, _ in pieces for x in (a, b)}
    cage = set(graph) - set().union(*(atoms for _, _, atoms in pieces))
    fusions = [(a, b, (rings[i].prefix, i)) for i, (a, b, _) in enumerate(pieces)]
    found = modified_candidate(mol, cage, fusions, fusion_atoms, _ranker(rings))
    return None if found is None else found.text(_renderer(rings))


def _pieces(mol, graph):
    """(atom, atom, atoms) for each group of atoms outside the cage that is bonded only to the two ends of one bond."""
    found = []
    for a, b in graph.edges:
        if graph.degree(a) < 3 or graph.degree(b) < 3 or any(mol.GetAtomWithIdx(x).GetAtomicNum() != 6 for x in (a, b)):
            continue
        parts = sorted(nx.connected_components(graph.subgraph(set(graph) - {a, b})), key=len)
        for part in parts[:-1]:
            if len(part) <= _MAX_PIECE and part & set(graph[a]) and part & set(graph[b]):
                found.append((a, b, frozenset(part)))
    pieces = [
        p for p in found
        if not any(q[2] | {q[0], q[1]} > p[2] | {p[0], p[1]} for q in found)
    ]
    used = set()
    for a, b, atoms in pieces:
        members = atoms | {a, b}
        if used & members or any(
            len(set(graph[x]) & atoms) != 1 or mol.GetAtomWithIdx(x).GetTotalNumHs() or graph.degree(x) not in (3, 4)
            for x in (a, b)
        ):
            return []
        used |= members
    return pieces


def _ring(mol, kekule, graph, a, b, atoms):
    """The fused component with the numberings that give its fusion atoms and its indicated hydrogen the lowest locants."""
    members = sorted(atoms | {a, b})
    index = {x: i for i, x in enumerate(members)}
    sub = skeleton(
        [mol.GetAtomWithIdx(x).GetSymbol() for x in members],
        [(index[x], index[y]) for x, y in graph.edges if x in index and y in index],
    )
    ring_info = sub.GetRingInfo()
    if any(not ring_info.NumAtomRings(i) for i in range(len(members))):
        return None
    if any(kekule.GetBondBetweenAtoms(x, n).GetBondType() == Chem.BondType.DOUBLE for x in (a, b) for n in graph[x] if n in atoms):
        return None
    try:
        components = identify(sub)
    except UnsupportedStructure:
        return None
    if not components or len({c.bracketed_prefix() for c in components}) > 1:
        return None
    component = components[0]

    def saturated(x):
        atom = mol.GetAtomWithIdx(x)
        return (
            atom.GetSymbol() not in _DONORS
            and atom.GetTotalNumHs() > 0
            and all(kekule.GetBondBetweenAtoms(x, n).GetBondType() != Chem.BondType.DOUBLE for n in graph[x] if n in index)
        )

    options = []
    for numbering in component.numberings:
        locant = {x: numbering[index[x]] for x in members}
        positions = sorted((locant[x] for x in atoms if saturated(x)), key=_locant_key)
        odd = len(positions) % 2
        indicated, hydro = tuple(positions[:odd]), tuple(positions[odd:])
        key = (
            tuple(sorted(_locant_key(locant[x]) for x in (a, b))),
            tuple(_locant_key(p) for p in indicated),
            tuple(_locant_key(p) for p in hydro),
        )
        options.append((key, (locant[a], locant[b], indicated, hydro)))
    best = min(key for key, _ in options)
    kind = component.kind in ("mono", "attached_only")
    return _Ring(
        prefix=component.prefix,
        bracketed=component.bracketed_prefix(),
        cite=not (kind and all(e == "C" for e in component.elements)),
        compound=bool(component.hetero_text)
        or (kind and component.name not in _RETAINED_PREFIX and any(e != "C" for e in component.elements)),
        options=[option for key, option in options if key == best],
    )


def _placed(rings, numbering, edit):
    """The fused rings grouped by prefix in alphanumerical order, each group by its lower fullerene locant; every ring
    with the locants of its fusion atoms cited in the direction of the fullerene locants."""
    groups = {}
    for atom_a, atom_b, (_, i) in edit.fusions:
        ring = rings[i]
        ascending = numbering[atom_a] < numbering[atom_b]
        cage = tuple(sorted((numbering[atom_a], numbering[atom_b])))

        def cited(option):
            first, second = (option[0], option[1]) if ascending else (option[1], option[0])
            return _locant_key(first), _locant_key(second)

        option = min(ring.options, key=cited)
        locants = (option[0], option[1]) if ascending else (option[1], option[0])
        groups.setdefault((ring.prefix, ring.bracketed), []).append((ring, cage, locants, option))
    return [sorted(groups[key], key=lambda item: item[1]) for key in sorted(groups)]


def _ranker(rings):
    def rank(numbering, edit):
        return tuple(_locant_key(l) for group in _placed(rings, numbering, edit) for item in group for l in item[2])

    return rank


def _renderer(rings):
    def render(numbering, edit):
        groups = _placed(rings, numbering, edit)
        prime, texts, indicated, hydro = 0, [], [], []
        for members in groups:
            sets = []
            for ring, cage, ring_locants, option in members:
                prime += 1
                mark = _PRIMES[prime] if prime < len(_PRIMES) else "′" * prime
                if ring.cite:
                    sets.append(f"{ring_locants[0]}{mark},{ring_locants[1]}{mark}:{cage[0]},{cage[1]}")
                else:
                    sets.append(f"{cage[0]},{cage[1]}")
                indicated += [((_locant_key(p)[0], 0, prime), f"{p}{mark}H") for p in option[2]]
                hydro += [((_locant_key(p), prime), f"{p}{mark}") for p in option[3]]
            ring = members[0][0]
            name = ring.bracketed
            if len(members) > 1:
                name = multiplying_prefix(len(members), compound=True) + f"({name})" if ring.compound else multiplying_prefix(len(members)) + name
            texts.append(f"{name}[{(';' if ring.cite else ':').join(sets)}]")
        hydro_text = ""
        if hydro:
            locants = [text for _, text in sorted(hydro)]
            hydro_text = ",".join(locants) + "-" + (multiplying_prefix(len(locants)) if len(locants) > 1 else "") + "hydro"
        return hydro_text, indicated, "".join(texts)

    return render
