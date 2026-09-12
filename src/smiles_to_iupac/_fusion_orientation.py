"""P-25.3.2.3.3 criterion (a) -- "maximum number of rings in a horizontal
row" -- generalized to any ortho-fused, non-peri-fused all-carbon aromatic
ring system whose ring-fusion graph is a tree (catacondensed chains, plus
branched systems like triphenylene) (https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf).

Model: each ring-fusion bond has a direction, in units of 60 degrees,
defined only up to one global rotation and modulo 6 (opposite directions
are the same lattice axis, `direction % 3`). Directions are assigned by
picking an arbitrary bond as direction 0, then propagating outward through
the ring-fusion tree: at each ring, every other incident bond's direction
is the already-known reference bond's direction plus a signed turn,
`(ref_edge - other_edge) % 6 - 3`, where `ref_edge`/`other_edge` are each
bond's position (0-5) on that ring's own hexagon, via `_edge_index`. This
generalizes the old chain-only "signed turn" model (diff 3 = 0 turn/
'straight', diff 4/2 = +-1 = 'bent') to a ring with more than 2 fusion
bonds: a third bond on adjacent hexagon edges (diff 1 or 5, e.g. a
branching center like triphenylene's) gets turn +-2, still consistent with
the same formula. A ring is "in" a candidate horizontal row (one of the 3
lattice axes) if any of its incident bond directions lies on that axis --
this is what lets a 'bent' ring still count as the last ring of a row it
bends away from, matching phenanthrene (row = 2, not 1), and lets a
branching ring's three neighbors each be checked independently.

Validated directly against the primary source's own worked examples
(P-25.3.2.3.3, P-25.3.2.4): anthracene = 3, phenanthrene = 2, tetraphene/
benzo[a]anthracene = 3, chrysene = 2 -- plus benzo[c]phenanthrene as a
check specifically for the signed-turn distinction (see the chain-only
version's history), and triphenylene = 2 as the first branching case,
matching `_triphenylene_fusion.py`'s own documented finding that chrysene
and triphenylene tie on every criterion up through (a) itself (the module
docstring there: "Chrysene and triphenylene tie on every criterion through
(f) ... down to (g)"), which is only possible if both compute the same row
count.

Scope, still narrow: the ring-fusion graph must be a tree (no cycle) --
peri-fused systems (e.g. pyrene, fluoranthene: an atom shared by three or
more rings, which also means the ring-fusion graph has a cycle rather than
just branching) remain out of scope and raise `UnsupportedStructure`; a
future step would need genuine 2D coordinate placement (not just a
direction per bond) to validate a peri-fused ring-fusion cycle closes
consistently. This module does not implement criteria (b)/(c)/(d) (the
quadrant-counting tie-breaks used when (a) alone doesn't decide -- needed
to actually resolve the chrysene vs triphenylene tie), and its result is
not yet wired into any naming path -- see
`_aromatic.py`/`_triphenylene_fusion.py`/`_chrysene_fusion.py`'s own
docstrings for where a real seniority decision still needs this.
"""

from ._common import UnsupportedStructure, adjacency, ring_cycle
from ._aromatic import (
    _atom_ring_membership,
    _edge_index,
    _ring_adjacency,
    find_aromatic_fused_core,
)


def _assign_bond_directions(graph, atom_rings, adj, fusion_bonds_by_pair, n):
    """Direction (mod 6, up to one global rotation) for every ring-fusion
    bond in a ring-fusion tree of `n` rings. Raises `UnsupportedStructure`
    if the ring-fusion graph is not a tree (a cycle means peri-fusion,
    e.g. pyrene)."""
    if n == 1:
        return {}
    total_edges = sum(len(v) for v in adj.values()) // 2
    if total_edges != n - 1:
        raise UnsupportedStructure(
            "a cyclic ring-fusion arrangement (e.g. a peri-fused system "
            "like pyrene or fluoranthene) is not supported yet (see "
            "P-25.3.2.3.3)"
        )

    direction = {}
    root = 0
    ref_neighbor = next(iter(adj[root]))
    direction[frozenset((root, ref_neighbor))] = 0
    processed_ref = {root: ref_neighbor, ref_neighbor: root}
    visited = {root, ref_neighbor}
    queue = [root, ref_neighbor]
    while queue:
        current = queue.pop(0)
        parent = processed_ref[current]
        cycle = ring_cycle(graph, list(atom_rings[current]))
        ref_edge = _edge_index(cycle, *fusion_bonds_by_pair[frozenset((current, parent))])
        for neighbor in adj[current]:
            if neighbor == parent:
                continue
            if neighbor in visited:
                raise UnsupportedStructure(
                    "a cyclic ring-fusion arrangement (e.g. a peri-fused "
                    "system like pyrene or fluoranthene) is not supported "
                    "yet (see P-25.3.2.3.3)"
                )
            visited.add(neighbor)
            bond = fusion_bonds_by_pair[frozenset((current, neighbor))]
            edge = _edge_index(cycle, *bond)
            diff = (ref_edge - edge) % 6
            direction[frozenset((current, neighbor))] = (
                direction[frozenset((current, parent))] + (diff - 3)
            )
            processed_ref[neighbor] = current
            queue.append(neighbor)

    if len(visited) != n:
        raise UnsupportedStructure(
            "a disconnected ring-fusion arrangement is not supported (see "
            "P-25.3.2.3.3)"
        )
    return direction


def _rings_on_axis(adj, direction, residue):
    count = 0
    for ring_idx, neighbors in adj.items():
        for neighbor in neighbors:
            if direction[frozenset((ring_idx, neighbor))] % 3 == residue:
                count += 1
                break
    return count


def max_rings_in_horizontal_row(adj, direction, n):
    """P-25.3.2.3.3 criterion (a): the largest number of rings that lie on
    one of the 3 possible horizontal-row axes, for a ring-fusion tree of
    `n` rings (`direction`: this tree's `_assign_bond_directions` result)."""
    if n <= 1:
        return n
    return max(_rings_on_axis(adj, direction, residue) for residue in range(3))


def count_rings_in_horizontal_row(mol) -> int:
    """P-25.3.2.3.3 criterion (a) for a plain all-carbon aromatic mancude
    ring system: the number of rings in the orientation's horizontal row.
    Raises `UnsupportedStructure` for a peri-fused or disconnected
    ring-fusion arrangement (see module docstring)."""
    core = find_aromatic_fused_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "not a plain all-carbon aromatic mancude ring system (see P-25.3.1.3)"
        )
    atom_rings, ring_atom_sets, fusion_bond_idxs = core
    n = len(atom_rings)
    membership = _atom_ring_membership(atom_rings)
    if any(len(rings) >= 3 for rings in membership.values()):
        raise UnsupportedStructure(
            "peri-fused aromatic ring systems (an atom shared by three or "
            "more rings) are not supported yet (see P-25.3.2.3.3)"
        )
    adj, fusion_bonds_by_pair = _ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol)
    graph = adjacency(mol)
    direction = _assign_bond_directions(graph, atom_rings, adj, fusion_bonds_by_pair, n)
    return max_rings_in_horizontal_row(adj, direction, n)
