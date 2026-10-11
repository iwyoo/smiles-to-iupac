"""Rings fused to a fullerene (P-27.6.1) and bridges between its atoms (P-27.6.2): the atoms outside the cage are groups
bonded to the rest through two atoms only; a group bonded to two adjacent cage atoms is a fused ring, any other group is a
bridge of the cage with its fused rings. The cage that remains is numbered as a modified fullerene."""

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure
from ._fullerene_bridged import Bridges
from ._fullerene_fused import layout, ring_component, ring_key, ring_positions, ring_text
from ._fullerene_modified import modified_candidate

_MAX_ATOMS = 160
_MAX_GROUP = 40
_MAX_ROUNDS = 6
_MIN_CAGE_CARBONS = 24


def name_fullerene_derivative(mol):
    """Name of `mol` as a modified fullerene, with the rings fused to it and the bridges it carries, else None. A cage
    that needs 'cyclo' or 'seco' only to take in the atoms of a ring or bridge is named with it as prefix instead."""
    plain = modified_candidate(mol)
    if plain is not None and not (plain.edit.seco or plain.edit.cyclo):
        return plain.text()
    try:
        attached = _attached_name(mol)
    except UnsupportedStructure:
        attached = None
    if attached is not None:
        return attached
    return None if plain is None else plain.text()


def _attached_name(mol):
    ring_info = mol.GetRingInfo()
    if mol.GetNumAtoms() > _MAX_ATOMS or ring_info.NumRings() < 11 or sum(len(r) == 5 for r in ring_info.AtomRings()) < 5:
        return None
    if sum(a.GetAtomicNum() == 6 and a.GetDegree() >= 3 and not a.GetTotalNumHs() for a in mol.GetAtoms()) < _MIN_CAGE_CARBONS:
        return None
    graph = nx.Graph((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds())
    tried = []
    for adjacent_only in (False, True):
        groups = _peel(graph, adjacent_only)
        if not groups or {g[2] for g in groups} in tried:
            continue
        tried.append({g[2] for g in groups})
        named = _compose(mol, graph, groups)
        if named is not None:
            return named
    return None


def _compose(mol, graph, groups):
    cage = set(graph) - set().union(*(atoms for _, _, atoms in groups))
    if len(cage) < _MIN_CAGE_CARBONS or not all(_attaches_once(mol, graph, x, groups) for x in cage if any(n not in cage for n in graph[x])):
        return None
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return None
    rings, bridges = [], []
    for u, v, atoms in groups:
        if u in cage and v in cage and graph.has_edge(u, v):
            ring = ring_component(mol, kekule, graph, u, v, atoms)
            if ring is None:
                return None
            rings.append(ring)
        else:
            bridges.append((u, v, atoms))
    fused = cage.union(*(r.atoms for r in rings))
    if any(u not in fused or v not in fused for u, v, _ in bridges):
        return None
    bridged = Bridges(mol, graph, fused, [atoms for _, _, atoms in bridges]) if bridges else None
    attach = {x for x in cage if any(n not in cage for n in graph[x])}
    fusions = [(r.a, r.b, (0, r.prefix, i)) for i, r in enumerate(rings)]
    fusions += [(u, v, (1, "", j)) for j, (u, v, _) in enumerate(bridges) if u in cage and v in cage]

    def placed_for(numbering, edit):
        groups_ = layout(rings, numbering, edit)
        position = {a: numbering[g] for a, g in edit.mapping.items()}
        position.update(ring_positions(groups_))
        position[-1] = max(numbering.values())  # bridge atoms are numbered on from the highest locant of the cage
        return groups_, position

    def rank(numbering, edit):
        groups_, position = placed_for(numbering, edit)
        bridge_key = bridged.assembled(position)[0][0] if bridged else ()
        return bridge_key, ring_key(groups_)

    def render(numbering, edit):
        groups_, position = placed_for(numbering, edit)
        hydro, indicated, fusion = ring_text(groups_)
        bridge = bridged.assembled(position)[0][1] if bridged else ""
        joiner = "-" if bridge and fusion and fusion[0] in "[0123456789" else ""
        return hydro, indicated, bridge + joiner + fusion

    found = modified_candidate(mol, cage, fusions, attach, rank)
    return None if found is None else found.text(render)


def _attaches_once(mol, graph, x, groups):
    """A cage atom that carries a group is a carbon without hydrogen of three or four bonds, bonded to one atom of it."""
    atom = mol.GetAtomWithIdx(x)
    outside = [n for n in graph[x] if any(n in atoms for _, _, atoms in groups)]
    return atom.GetAtomicNum() == 6 and not atom.GetTotalNumHs() and graph.degree(x) in (3, 4) and len(outside) == 1


def _peel(graph, adjacent_only):
    """(atom, atom, atoms) for each group of atoms bonded to the rest through two atoms only, taken off in rounds so that
    a group between rings goes before the rings; a group that contains another is taken whole. The boundary of a
    fullerene fragment looks like a bridge between two atoms, so `adjacent_only` keeps the groups of adjacent atoms."""
    work = graph.copy()
    found = []
    for _ in range(_MAX_ROUNDS):
        cuts = [c for c in _cuts(work) if not adjacent_only or work.has_edge(c[0], c[1])]
        if not cuts:
            return found
        closure = [c[2] | {c[0], c[1]} for c in cuts]
        taken = [c for c, outer in zip(cuts, closure) if not any(outer < other for other in closure)]
        members = [x for c in taken for x in c[2]]
        if len(members) != len(set(members)):
            return []
        found += taken
        work.remove_nodes_from(members)
    return []


def _cuts(graph):
    found = {}
    for u in graph:
        rest = graph.copy()
        rest.remove_node(u)
        for v in nx.articulation_points(rest):
            parts = sorted(nx.connected_components(rest.subgraph(set(rest) - {v})), key=len)
            for part in parts[:-1]:
                if len(part) <= _MAX_GROUP and part & set(graph[u]) and part & set(graph[v]):
                    found[(min(u, v), max(u, v), frozenset(part))] = True
    return list(found)
