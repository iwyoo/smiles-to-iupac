"""P-25.3.2.3.3 criteria (a) and (b) -- "maximum number of rings in a
horizontal row" and "maximum number of rings in the upper right quadrant"
-- generalized to any ortho-fused, non-peri-fused all-carbon aromatic ring
system whose ring-fusion graph is a tree (catacondensed chains, plus
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

Criterion (b) extends the same bond-direction data to actual 2D ring
coordinates (`_ring_positions`), scaled so every position -- and every
row's center point, even a bond midpoint -- lands on an integer (see
`_HEX_DELTA`), so every "is this ring exactly on an axis" check is exact
integer comparison, never a float tolerance. A ring exactly on one axis
counts as a half (P-25.3.2.3.3(b): "those rings that are divided by an
axis are considered as two halves"), and a ring at the intersection of
both axes as a quarter. Validated against phenanthrene = 1.5, the primary
source's own worked example -- and, once computed, chrysene and
triphenylene tie at 2.5 as well as at criterion (a)'s 2, one further data
point supporting `_triphenylene_fusion.py`'s existing guess that the two
tie all the way through P-25.3.2.4's seniority list and fall to plain
alphabetical order, rather than being decided by these quadrant criteria
at all -- worth re-checking once (c)/(d) exist too, since if the tie
really does persist through all four, chasing it further with (c)/(d) is
the wrong next step for resolving *this specific* comparison (though (c)/
(d) are still needed generally, for other ties and for a single
structure's own preferred-orientation choice).

Scope, still narrow: the ring-fusion graph must be a tree (no cycle) --
peri-fused systems (e.g. pyrene, fluoranthene: an atom shared by three or
more rings, which also means the ring-fusion graph has a cycle rather than
just branching) remain out of scope and raise `UnsupportedStructure`; a
future step would need genuine 2D coordinate placement for a *cyclic*
ring-fusion graph (this module's 2D placement so far only walks a tree) to
validate it closes consistently. This module does not implement criteria
(c)/(d) (the remaining quadrant-counting tie-breaks), and its result is
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


def _prepare(mol):
    """Shared preamble for every P-25.3.2.3.3 criterion: validate `mol` is
    a plain all-carbon aromatic mancude ring system with a tree-shaped
    (non-peri-fused) ring-fusion graph, and return (atom_rings, adj,
    fusion_bonds_by_pair, direction, n)."""
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
    return atom_rings, adj, fusion_bonds_by_pair, direction, n


def count_rings_in_horizontal_row(mol) -> int:
    """P-25.3.2.3.3 criterion (a) for a plain all-carbon aromatic mancude
    ring system: the number of rings in the orientation's horizontal row.
    Raises `UnsupportedStructure` for a peri-fused or disconnected
    ring-fusion arrangement (see module docstring)."""
    _atom_rings, adj, _fusion_bonds_by_pair, direction, n = _prepare(mol)
    return max_rings_in_horizontal_row(adj, direction, n)


# Hex-lattice unit vector for each bond direction (0-5), scaled by 2 along
# x and by 2/sqrt(3) along y so every ring center -- and every row's
# center point, even a bond midpoint -- lands on an *integer* coordinate.
# This keeps every "is this ring exactly on an axis" check in criterion
# (b)/(c)/(d) an exact integer comparison, never a float tolerance.
_HEX_DELTA = {0: (2, 0), 1: (1, 1), 2: (-1, 1), 3: (-2, 0), 4: (-1, -1), 5: (1, -1)}


def _ring_positions(adj, direction, winning_residue, n):
    """Integer (x, y) position for every ring (see `_HEX_DELTA`), oriented
    so bonds on `winning_residue`'s axis (criterion (a)'s winning axis)
    run horizontally. Positions are relative to an arbitrary ring (index
    0) at the origin, not yet shifted to the row's center -- see
    `_row_center`."""
    position = {0: (0, 0)}
    visited = {0}
    queue = [0]
    while queue:
        current = queue.pop(0)
        for neighbor in adj[current]:
            if neighbor in visited:
                continue
            visited.add(neighbor)
            raw = direction[frozenset((current, neighbor))]
            dx, dy = _HEX_DELTA[(raw - winning_residue) % 6]
            cx, cy = position[current]
            position[neighbor] = (cx + dx, cy + dy)
            queue.append(neighbor)
    return position


def _row_center(adj, direction, position, winning_residue, n):
    """P-25.3.2.3.3(b)'s "center of the horizontal row": the central
    common bond's midpoint if the row has an even number of rings, the
    central ring's own position if odd. The row is the longest
    connected run of rings joined only by bonds on `winning_residue`'s
    axis (every such bond keeps y fixed and only changes x by +-2, so the
    row's members share one y and are ordered by x)."""
    axis_adj = {i: set() for i in range(n)}
    for ring_idx, neighbors in adj.items():
        for neighbor in neighbors:
            if direction[frozenset((ring_idx, neighbor))] % 3 == winning_residue:
                axis_adj[ring_idx].add(neighbor)

    seen = set()
    best_component = []
    for start in range(n):
        if start in seen:
            continue
        component = []
        stack = [start]
        seen.add(start)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in axis_adj[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        if len(component) > len(best_component):
            best_component = component

    row = sorted(best_component, key=lambda i: position[i][0])
    m = len(row)
    if m % 2 == 0:
        left, right = position[row[m // 2 - 1]], position[row[m // 2]]
        return ((left[0] + right[0]) // 2, left[1])
    center_ring = position[row[m // 2]]
    return center_ring


def rings_in_upper_right_quadrant(mol) -> float:
    """P-25.3.2.3.3 criterion (b): the number of rings in the upper right
    quadrant relative to the horizontal row's center, per criterion (a).
    A ring exactly on one axis counts as a half (attributed to whichever
    of the two quadrants that axis's side actually borders); a ring at
    the exact intersection of both axes counts as a quarter. Raises
    `UnsupportedStructure` under the same conditions as
    `count_rings_in_horizontal_row`.

    Validated against the primary source's own worked example
    (P-25.3.2.3.3(b)): phenanthrene = 1.5 (phenalene, the other half of
    that example, is peri-fused and out of scope here)."""
    _atom_rings, adj, _fusion_bonds_by_pair, direction, n = _prepare(mol)
    if n <= 1:
        return 0.25 if n == 1 else 0.0
    winning_residue = max(range(3), key=lambda r: _rings_on_axis(adj, direction, r))
    position = _ring_positions(adj, direction, winning_residue, n)
    center_x, center_y = _row_center(adj, direction, position, winning_residue, n)

    total = 0.0
    for ring_idx in range(n):
        x = position[ring_idx][0] - center_x
        y = position[ring_idx][1] - center_y
        if x == 0 and y == 0:
            total += 0.25
        elif x == 0:
            if y > 0:
                total += 0.5
        elif y == 0:
            if x > 0:
                total += 0.5
        elif x > 0 and y > 0:
            total += 1.0
    return total
