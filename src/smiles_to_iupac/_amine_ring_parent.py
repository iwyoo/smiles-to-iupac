"""Polyamines whose parent is a benzene ring (P-62.2.4.1.3, P-44.1.1, P-45.2.1): the ring bearing the most amine
nitrogens is the parent, every other group on those nitrogens is an N-locanted prefix, other rings and chains are
cited as substituent prefixes. N1-(4-aminophenyl)-N4-phenylbenzene-1,4-diamine."""

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    lowest_locant_set,
    multiplied_word,
    ring_cycle,
    separate_aromatic_monocycles,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_ring


def chain_amine_count(mol, graph):
    """Most amine nitrogens one acyclic carbon chain can carry as suffixes."""
    chain_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
    nitrogens = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    best = 0

    def walk(path, seen):
        nonlocal best
        on_path = set(path)
        carried = sum(1 for n in nitrogens if sum(c in on_path for c in graph[n]) == 1)
        best = max(best, carried)
        for nxt in graph[path[-1]]:
            if nxt in chain_atoms and nxt not in seen:
                walk(path + [nxt], seen | {nxt})

    for start in chain_atoms:
        walk([start], {start})
    return best


def _plain_rings(mol, graph):
    rings = separate_aromatic_monocycles(mol, graph)
    if rings is None:
        ring_atoms = [set(r) for r in mol.GetRingInfo().AtomRings()]
        rings = ring_atoms if len(ring_atoms) == 1 else None
    if rings is None or not all(is_plain_benzene_ring(mol, r) for r in rings):
        raise UnsupportedStructure("rings other than separate benzene rings alongside several amines are not supported yet")
    return rings


def _check_amine_nitrogens(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 7, 9, 17, 35, 53):
            raise UnsupportedStructure("heteroatoms other than amine nitrogen and halogens are not supported here")
        if atom.GetFormalCharge() or atom.GetIsotope() or atom.GetNumRadicalElectrons():
            raise UnsupportedStructure("charged, radical or isotopically modified atoms are not supported here")
        if atom.GetAtomicNum() == 7 and (
            atom.GetIsAromatic()
            or any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors())
            or any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds())
        ):
            raise UnsupportedStructure("a nitrogen that is not a plain amine nitrogen is not supported here")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0 and not bond.GetIsAromatic():
            raise UnsupportedStructure("non-aromatic multiple bonds alongside several amines are not supported here")


def _parent_nitrogens(graph, ring_atoms, nitrogens):
    return [n for n in nitrogens if sum(c in ring_atoms for c in graph[n]) == 1]


def _candidate(mol, graph, halogens, ring_order, parent_nitrogens):
    from ._amine import _add_n_names

    n_names_by_nitrogen = {}
    for n in parent_nitrogens:
        n_names_by_nitrogen[n] = [
            name_branch(graph, c, n, halogens, mol=mol) for c in graph[n] if c not in ring_order
        ]
    count = len(parent_nitrogens)
    best = None
    for start in range(len(ring_order)):
        rotated = ring_order[start:] + ring_order[:start]
        for numbering in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(numbering)}
            locant_of = {n: position_of[next(c for c in graph[n] if c in position_of)] for n in parent_nitrogens}
            amine_locants = sorted(locant_of.values())
            ring_substituents = substituents_for_ring(graph, numbering, halogens, set(parent_nitrogens), mol=mol)
            grouped = group_substituents(ring_substituents)
            n_names, n_locants = [], []
            for n in sorted(parent_nitrogens, key=lambda n: locant_of[n]):
                for entry in sorted(n_names_by_nitrogen[n]):
                    n_names.append(entry)
                    n_locants.append(None if count == 1 else locant_of[n])
            locant_set, total, citation = substituent_locant_set_and_citation(grouped)
            prefix = format_substituent_prefixes(_add_n_names(grouped, n_names, n_locants))
            if count == 1:
                name = f"{prefix}aniline"
            else:
                name = f"{prefix}benzene-{','.join(map(str, amine_locants))}-{multiplied_word(count, 'amine')}"
            total_prefixes = total + len(n_names)
            key = (lowest_locant_set(amine_locants), locant_set, citation, name)
            if best is None or key < best[0]:
                best = (key, name, total_prefixes)
    return best


def name_ring_parent_polyamine(mol):
    if specified_stereocenters(mol):
        raise UnsupportedStructure("a specified stereocenter alongside a polyamine is not supported yet")
    _check_amine_nitrogens(mol)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    nitrogens = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    rings = _plain_rings(mol, graph)
    options = []
    for ring in rings:
        carried = _parent_nitrogens(graph, ring, nitrogens)
        if carried:
            options.append((ring, carried))
    if not options:
        return None
    most = max(len(carried) for _, carried in options)
    if chain_amine_count(mol, graph) > most:
        return None
    best = None
    for ring, carried in options:
        if len(carried) != most:
            continue
        found = _candidate(mol, graph, halogens, ring_cycle(graph, list(ring)), carried)
        key = (-found[2], found[0])
        if best is None or key < best[0]:
            best = (key, found[1])
    return best[1]
