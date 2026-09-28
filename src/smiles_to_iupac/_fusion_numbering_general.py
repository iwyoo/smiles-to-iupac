"""Atom-level P-25.3.3.1 peripheral numbering for a tree of ortho-fused
mancude rings, built on `_fusion_orientation_general.py`'s ring-level
starting-ring engine, plus P-25.3.3.1.2's heteroatom/indicated-H
tie-break. Exactly reproduces the already-verified all-hexagon
benzo[g]quinoline numbering. NOT yet trustworthy for a tree containing a
non-hexagonal ring: cross-checking against the two other real numberings
this session shipped by hand (1H-/3H-cyclopenta[a]naphthalene, 1H-
cyclopenta[b]naphthalene) shows the ring-level engine picks the wrong
starting ring for them -- its Fraction-of-a-turn model treats every ring
size as a regular polygon, but P-25.3.2.3.1's actual permitted shapes are
fixed, non-regular drawings (see issue filed against this finding).
Deliberately NOT wired into `core.py`: it would silently regress those
two already-correct real cases.
"""

from collections import defaultdict

from ._fusion_orientation_general import assign_bond_directions_general, starting_ring_general


def _ring_bond_cycle(mol, atoms):
    adjacency = defaultdict(list)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in atoms and b in atoms:
            adjacency[a].append(b)
            adjacency[b].append(a)
    start = min(atoms)
    cycle = [start]
    prev, cur = None, start
    while True:
        nxt = next(x for x in adjacency[cur] if x != prev)
        if nxt == start:
            break
        cycle.append(nxt)
        prev, cur = cur, nxt
    return [frozenset((cycle[i], cycle[(i + 1) % len(cycle)])) for i in range(len(cycle))]


def _ring_graph(mol):
    """{ring_idx: atom set}, adj, edge_index_of, ring_sizes for
    `_fusion_orientation_general.py`'s interface -- or None if the fusion
    graph isn't a simple ortho-fused tree (peri-fusion, spiro, etc.)."""
    ring_info = mol.GetRingInfo()
    atom_rings = [set(r) for r in ring_info.AtomRings()]
    n = len(atom_rings)
    bond_cycles = [_ring_bond_cycle(mol, atoms) for atoms in atom_rings]
    adj = {i: set() for i in range(n)}
    edge_index_of = {}
    for i in range(n):
        for j in range(i + 1, n):
            shared = atom_rings[i] & atom_rings[j]
            if len(shared) != 2:
                continue
            bond = frozenset(shared)
            if bond not in bond_cycles[i] or bond not in bond_cycles[j]:
                continue
            adj[i].add(j)
            adj[j].add(i)
            edge_index_of[(i, j)] = bond_cycles[i].index(bond)
            edge_index_of[(j, i)] = bond_cycles[j].index(bond)
    ring_sizes = [len(atoms) for atoms in atom_rings]
    return atom_rings, adj, edge_index_of, ring_sizes, n


def _periphery_cycle_and_fusion_atoms(mol):
    ring_info = mol.GetRingInfo()
    bond_count = defaultdict(int)
    for bond_ring in ring_info.BondRings():
        for bond_idx in bond_ring:
            bond_count[bond_idx] += 1
    interior_bonds = {b for b, c in bond_count.items() if c >= 2}
    fusion_atoms = set()
    for bond_idx in interior_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        fusion_atoms.add(bond.GetBeginAtomIdx())
        fusion_atoms.add(bond.GetEndAtomIdx())

    adjacency = defaultdict(list)
    for bond in mol.GetBonds():
        if bond.GetIdx() in interior_bonds:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adjacency[a].append(b)
        adjacency[b].append(a)
    start = next(iter(adjacency))
    cycle = [start]
    prev, cur = None, start
    while True:
        nxt = next(x for x in adjacency[cur] if x != prev)
        if nxt == start:
            break
        cycle.append(nxt)
        prev, cur = cur, nxt
    return cycle, fusion_atoms


def _numbering_from(cycle, fusion_atoms, start_i, direction):
    n = len(cycle)
    locants = {}
    counter = 0
    letter_ord = 0
    prev_was_fusion = False
    for step in range(n):
        atom = cycle[(start_i + direction * step) % n]
        if atom in fusion_atoms:
            letter_ord = letter_ord + 1 if prev_was_fusion else 1
            locants[atom] = f"{counter}{chr(ord('a') + letter_ord - 1)}"
            prev_was_fusion = True
        else:
            counter += 1
            locants[atom] = str(counter)
            prev_was_fusion = False
    return locants


def general_peripheral_numbering(mol):
    """{atom_idx: locant_str} for a tree of ortho-fused mancude rings of
    any size, or None if the fusion graph isn't a simple ortho-fused tree.
    Heteroatom-lowest-locants is the only P-25.3.3.1.2 tie-break applied;
    ties beyond that resolve by ring/atom traversal order, not further
    Blue Book sub-criteria (P-25.3.3.1.2(d)-(f) are not implemented)."""
    atom_rings, adj, edge_index_of, ring_sizes, n = _ring_graph(mol)
    if n < 2 or sum(len(v) for v in adj.values()) // 2 != n - 1:
        return None

    direction = assign_bond_directions_general(adj, edge_index_of, ring_sizes, n)
    start_rings = starting_ring_general(adj, direction, n)

    cycle, fusion_atoms = _periphery_cycle_and_fusion_atoms(mol)
    idx_of = {a: i for i, a in enumerate(cycle)}
    hetero = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 6]

    candidates = []
    for ring_idx in start_rings:
        for f_atom in atom_rings[ring_idx] & fusion_atoms:
            i = idx_of[f_atom]
            for step in (1, -1):
                neighbor = cycle[(i + step) % len(cycle)]
                if neighbor in atom_rings[ring_idx] and neighbor not in fusion_atoms:
                    candidates.append(_numbering_from(cycle, fusion_atoms, (i + step) % len(cycle), step))

    if not candidates:
        return None

    indicated_h = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 2]

    def sort_key(locants):
        hetero_key = sorted(int(locants[h].rstrip("abcdefgh")) for h in hetero)
        indicated_key = sorted(int(locants[h].rstrip("abcdefgh")) for h in indicated_h)
        return (hetero_key, indicated_key)

    return min(candidates, key=sort_key)
