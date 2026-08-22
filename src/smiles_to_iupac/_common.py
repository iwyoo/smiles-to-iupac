"""Shared graph utilities and validation for saturated hydrocarbon parent
hydrides, used by both the acyclic (`_acyclic.py`) and monocyclic (`_cyclic.py`)
naming modules.
"""

from rdkit import Chem


class UnsupportedStructure(NotImplementedError):
    pass


def validate_atoms_and_bonds(mol):
    """Structure-independent checks shared by every parent hydride kind: an
    all-carbon, single-fragment skeleton. Bond order (all single, or exactly one
    double/triple bond) and ring shape (none, one simple ring, or more) are
    checked separately by each naming module, since what's allowed there
    differs (see `_acyclic.py`, `_cyclic.py`, `_unsaturated.py`)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms are not supported yet (see P-21.2.3, skeletal "
                "replacement nomenclature)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )


def non_single_bonds(mol):
    """List of (begin_atom_idx, end_atom_idx, bond_order) for every bond whose
    order isn't 1.0 (single). Used to classify a molecule's degree of
    unsaturation for dispatch (see `core.py`)."""
    return [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0
    ]


def adjacency(mol):
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        graph[a].append(b)
        graph[b].append(a)
    return graph


def bfs(graph, start):
    dist = {start: 0}
    parent = {start: None}
    queue = [start]
    while queue:
        next_queue = []
        for node in queue:
            for neighbor in graph[node]:
                if neighbor not in dist:
                    dist[neighbor] = dist[node] + 1
                    parent[neighbor] = node
                    next_queue.append(neighbor)
        queue = next_queue
    return dist, parent


def path_between(parent, start, end):
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    return list(reversed(path))


def linear_branch(graph, root, coming_from):
    """Walk a branch outward; return its atom count, or None if it forks
    (a "compound" substituent, P-29.4, not yet supported)."""
    length = 1
    previous, current = coming_from, root
    while True:
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) == 0:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]
        length += 1


def lowest_locant_set(locants):
    return tuple(sorted(locants))
