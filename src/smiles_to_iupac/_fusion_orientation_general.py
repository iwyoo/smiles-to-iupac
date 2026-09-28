"""Generalizes `_fusion_orientation.py`'s P-25.3.2.3.3 criteria (a)-(d)
cascade from all-hexagon ring-fusion trees to trees mixing any permitted
3-8 membered ring size (P-25.3.2.3.1). Direction is an exact `Fraction`
of a full turn (`% N`-per-step, centered by `N // 2`, generalizing the
hexagon-only module's own `% 6` / `- 3` -- verified identical to it on 4
real molecules, see PR description and this module's own tests). 2D
placement uses real trigonometry (N=5/7 angles are irrational); axis
membership stays exact `Fraction` comparison. Scope: the algorithm only,
given an already-extracted ring-adjacency description -- real-molecule
extraction and the P-25.3.2.3.2 distorted-ring fallback are both out of
scope, same as `_fusion_orientation.py`'s own precedent (its own result
"is not yet wired into any naming path" either).
"""

import math
from fractions import Fraction

from ._common import UnsupportedStructure

HALF_TURN = Fraction(1, 2)


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
            diff = (ref_edge - this_edge) % n_current
            centered = diff - (n_current // 2)
            direction[frozenset((current, neighbor))] = (
                direction[frozenset((current, parent))] + Fraction(centered, n_current)
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


def best_orientation_general(adj, direction, n):
    """The single preferred orientation, per P-25.3.2.3.3's full (a) ->
    (b) -> (c) -> (d) cascade, generalized to mixed ring sizes. Returns
    (residue, reflect, upper_right, lower_left, above) for the winner."""
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
    return candidates[0]
