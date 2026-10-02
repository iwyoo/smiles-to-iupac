"""Generalizes `_fusion_orientation.py`'s P-25.3.2.3.3 orientation
criteria plus P-25.3.3.1.1's numbering-start rule to any mix of permitted
3-8 membered rings. `_PENTAGON_TABLE` encodes P-25.3.2.3.1's real (non-
regular) house shape, pixel-measured off `tmp/bluebook/p25_shapes/*.png`.
Verified against 2 of 3 real cyclopenta-naphthalene isomers.
"""

import math
from fractions import Fraction

from ._common import UnsupportedStructure

HALF_TURN = Fraction(1, 2)

# Index i = outward-normal angle of edge i for P-25.3.2.3.1's apex-up
# "house" pentagon (edge 0 = right vertical, then bottom/left-vertical/
# left-roof/right-roof), pixel-measured off the Blue Book's own diagram.
_PENTAGON_TABLE = [Fraction(0), Fraction(3, 4), Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)]
_SHAPE_TABLES = {5: _PENTAGON_TABLE}


def _relative_turn(ref_edge, this_edge, n_current):
    """Real diagram-derived edge geometry where a shape table exists
    (currently n=5 only); the old regular-polygon formula otherwise,
    unchanged (exact for the already-regular sizes 3/4/6)."""
    table = _SHAPE_TABLES.get(n_current)
    if table is not None:
        return (table[ref_edge] - table[this_edge] - HALF_TURN) % 1
    diff = (ref_edge - this_edge) % n_current
    centered = diff - (n_current // 2)
    return Fraction(centered, n_current)


def assign_bond_directions_general(adj, edge_index_of, ring_sizes, n):
    """Direction (exact `Fraction` of a full turn, up to one global
    rotation) for every ring-fusion bond in a tree of `n` rings, each of
    any permitted P-25.3.2.3.1 size. `adj`: {ring_idx: set(neighbor_idx)}.
    `edge_index_of`: {(ring_idx, neighbor_idx): int}, the fusion bond's
    position (0..ring_sizes[ring_idx]-1) in `ring_idx`'s own cyclic bond
    order, directed (may differ for the same bond as seen from each of
    its two rings, since they generally have different sizes/orderings).
    Raises `UnsupportedStructure` if the ring-fusion graph is not a tree
    (a cycle means peri-fusion)."""
    if n == 1:
        return {}
    total_edges = sum(len(v) for v in adj.values()) // 2
    if total_edges != n - 1:
        raise UnsupportedStructure(
            "a cyclic ring-fusion arrangement (peri-fusion) is not "
            "supported by this generalization yet (see P-25.3.2.3.3)"
        )

    direction = {}
    root = 0
    ref_neighbor = next(iter(adj[root]))
    direction[frozenset((root, ref_neighbor))] = Fraction(0)
    processed_ref = {root: ref_neighbor, ref_neighbor: root}
    visited = {root, ref_neighbor}
    queue = [root, ref_neighbor]
    while queue:
        current = queue.pop(0)
        parent = processed_ref[current]
        n_current = ring_sizes[current]
        ref_edge = edge_index_of[(current, parent)]
        for neighbor in adj[current]:
            if neighbor == parent:
                continue
            if neighbor in visited:
                raise UnsupportedStructure(
                    "a cyclic ring-fusion arrangement (peri-fusion) is not "
                    "supported by this generalization yet (see P-25.3.2.3.3)"
                )
            visited.add(neighbor)
            this_edge = edge_index_of[(current, neighbor)]
            direction[frozenset((current, neighbor))] = (
                direction[frozenset((current, parent))] + _relative_turn(ref_edge, this_edge, n_current)
            ) % 1
            processed_ref[neighbor] = current
            queue.append(neighbor)

    if len(visited) != n:
        raise UnsupportedStructure(
            "a disconnected ring-fusion arrangement is not supported (see "
            "P-25.3.2.3.3)"
        )
    return direction


def _rings_on_axis_general(adj, direction, residue):
    count = 0
    for ring_idx, neighbors in adj.items():
        for neighbor in neighbors:
            if direction[frozenset((ring_idx, neighbor))] % HALF_TURN == residue:
                count += 1
                break
    return count


def _axis_candidates(direction):
    """Every distinct axis class (`Fraction` in [0, 1/2)) that at least
    one fusion bond actually lies on -- generalizes the hexagon-only
    version's fixed 3-candidate search (residues 0,1,2 out of the 6
    possible hexagon directions) to however many distinct axis classes
    actually occur in this specific mixed-ring-size tree."""
    return sorted({d % HALF_TURN for d in direction.values()})


def max_rings_in_horizontal_row_general(adj, direction, n):
    """P-25.3.2.3.3 criterion (a), generalized: the largest number of
    rings that lie on any axis actually realized by a fusion bond in this
    tree of `n` rings."""
    if n <= 1:
        return n
    return max(_rings_on_axis_general(adj, direction, r) for r in _axis_candidates(direction))


_REFLECTIONS = ((1, 1), (-1, 1), (1, -1), (-1, -1))
_EPS = 1e-9


def _ring_positions_general(adj, direction, winning_residue, reflect=(1, 1)):
    """Float (x, y) position for every ring, oriented so bonds on
    `winning_residue`'s axis run horizontally. Real trigonometry (not
    exact integers) is required here -- see module docstring -- since
    N=5/7 edge angles are irrational multiples of pi. Every ring-to-ring
    hop is treated as unit length regardless of the two rings' actual
    sizes: criteria (a)-(d) are about topological row/quadrant placement,
    not metric ring geometry, and this reduces to the hexagon-only
    module's own (differently-scaled but proportionally equivalent)
    placement whenever every ring is N=6."""
    sx, sy = reflect
    position = {0: (0.0, 0.0)}
    visited = {0}
    queue = [0]
    while queue:
        current = queue.pop(0)
        for neighbor in adj[current]:
            if neighbor in visited:
                continue
            visited.add(neighbor)
            raw = direction[frozenset((current, neighbor))]
            theta = 2 * math.pi * float(raw - winning_residue)
            dx, dy = math.cos(theta), math.sin(theta)
            cx, cy = position[current]
            position[neighbor] = (cx + sx * dx, cy + sy * dy)
            queue.append(neighbor)
    return position


def _row_center_general(adj, direction, position, winning_residue, n):
    axis_adj = {i: set() for i in range(n)}
    for ring_idx, neighbors in adj.items():
        for neighbor in neighbors:
            if direction[frozenset((ring_idx, neighbor))] % HALF_TURN == winning_residue:
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
        return ((left[0] + right[0]) / 2, left[1])
    return position[row[m // 2]]


def _quadrant_counts_general(adj, direction, residue, reflect, n):
    position = _ring_positions_general(adj, direction, residue, reflect)
    center_x, center_y = _row_center_general(adj, direction, position, residue, n)
    upper_right = lower_left = above = 0.0
    for ring_idx in range(n):
        x = position[ring_idx][0] - center_x
        y = position[ring_idx][1] - center_y
        on_x_axis = abs(x) < _EPS
        on_y_axis = abs(y) < _EPS
        if on_x_axis and on_y_axis:
            upper_right += 0.25
            lower_left += 0.25
            above += 0.5
        elif on_x_axis:
            if y > 0:
                upper_right += 0.5
            else:
                lower_left += 0.5
        elif on_y_axis:
            if x > 0:
                upper_right += 0.5
            else:
                lower_left += 0.5
            above += 0.5
        else:
            if x > 0 and y > 0:
                upper_right += 1.0
            elif x < 0 and y < 0:
                lower_left += 1.0
            if y > 0:
                above += 1.0
    return upper_right, lower_left, above


def _all_best_orientations(adj, direction, n):
    """Every (residue, reflect, upper_right, lower_left, above) tied
    through the full (a)->(b)->(c)->(d) cascade -- an irregular ring size
    can make two tied candidates pick different rings as "uppermost",
    since criteria (a)-(d) never compare ring shapes."""
    residues = _axis_candidates(direction)
    best_a = max(_rings_on_axis_general(adj, direction, r) for r in residues)
    residues = [r for r in residues if _rings_on_axis_general(adj, direction, r) == best_a]
    candidates = [
        (residue, reflect) + _quadrant_counts_general(adj, direction, residue, reflect, n)
        for residue in residues
        for reflect in _REFLECTIONS
    ]
    max_b = max(c[2] for c in candidates)
    candidates = [c for c in candidates if abs(c[2] - max_b) < _EPS]
    min_c = min(c[3] for c in candidates)
    candidates = [c for c in candidates if abs(c[3] - min_c) < _EPS]
    max_d = max(c[4] for c in candidates)
    candidates = [c for c in candidates if abs(c[4] - max_d) < _EPS]
    return candidates


def best_orientation_general(adj, direction, n):
    """The single preferred orientation, per P-25.3.2.3.3's full (a) ->
    (b) -> (c) -> (d) cascade, generalized to mixed ring sizes. Returns
    (residue, reflect, upper_right, lower_left, above) for the winner --
    the first of possibly several criteria-tied candidates, see
    `_all_best_orientations`."""
    return _all_best_orientations(adj, direction, n)[0]


def starting_ring_general(adj, direction, n):
    """P-25.3.3.1.1's "uppermost, farthest right ring" rule, unioned
    across every criteria-(a)-(d)-tied orientation (see
    `_all_best_orientations`). Caller applies its own P-25.3.3.1.2
    tie-break among the returned candidates."""
    if n == 1:
        return [0]
    winners = set()
    for residue, reflect, _ur, _ll, _above in _all_best_orientations(adj, direction, n):
        position = _ring_positions_general(adj, direction, residue, reflect)
        max_y = max(y for _x, y in position.values())
        top_rings = [i for i in range(n) if abs(position[i][1] - max_y) < _EPS]
        max_x = max(position[i][0] for i in top_rings)
        winners.update(i for i in top_rings if abs(position[i][0] - max_x) < _EPS)
    return sorted(winners)
