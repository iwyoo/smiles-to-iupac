"""P-25.3.2.3.3 criterion (a) -- "maximum number of rings in a horizontal
row" -- generalized to catacondensed, non-branching ortho-fused aromatic
ring chains of arbitrary length (https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf).

`_aromatic.py`'s `_classify_shape` already tells each internal ring apart
as 'straight' or 'bent', but collapses the *direction* of a bend (two
possible geometries, `_edge_index` diff 2 vs 4) into that one label. That
loses information a general row-count needs: a chain that bends the same
way twice (e.g. chrysene, PubChem CID 9171, diffs [4, 4]) ends up with a
different horizontal-row count than one that bends, then bends back
(e.g. benzo[c]phenanthrene / [4]helicene, CID 9136, diffs [4, 2]), even
though both collapse to the identical ('bent', 'bent') shape tuple. So
this module recomputes the per-junction turn itself, signed (+1/-1 units
of 60 degrees; 0 for 'straight'), rather than consuming `_classify_shape`'s
already-collapsed output.

Model: each ring-fusion bond has a direction, in units of 60 degrees,
accumulated from each internal ring's own signed turn (diff==3 is 0 turn/
'straight'; diff==4 is +1 and diff==2 is -1, both 'bent' -- this sign is
internally consistent along one chain because every ring's own two fusion
bonds are compared via the same `ring_cycle` traversal direction, but is
not claimed to correspond to any fixed real-world "clockwise"/
"counterclockwise" outside this module). A ring is "in" a candidate
horizontal row (one of the 3 possible lattice axes, since a hexagon has 3
pairs of opposite edges) if either the bond direction entering it or the
bond direction leaving it lies on that axis (direction mod 3 == the
axis's residue) -- this is what lets a 'bent' ring still count as the
last ring of a row it bends away from, matching phenanthrene (row = 2,
not 1). Criterion (a)'s answer is the maximum count over the 3 axis
choices.

Validated directly against the primary source's own worked examples
(P-25.3.2.3.3): anthracene = 3, phenanthrene = 2, tetraphene/
benzo[a]anthracene = 3, chrysene = 2 -- plus benzo[c]phenanthrene as a
fifth check specifically for the signed-turn distinction above.

Scope, deliberately narrow: only catacondensed, non-branching,
non-peri-fused all-carbon aromatic chains (the same shape
`_ring_path_order` already restricts `_aromatic.py` to) -- branching
(e.g. triphenylene) and peri-fusion (e.g. pyrene, fluoranthene) are out of
scope for this module and raise `UnsupportedStructure`, same as
`_aromatic.py` does. This module does not implement criteria (b)/(c)/(d)
(the quadrant-counting tie-breaks used when (a) alone doesn't decide), and
its result is not yet wired into any naming path -- see
`_aromatic.py`/`_triphenylene_fusion.py`/`_chrysene_fusion.py`'s own
docstrings for where a real seniority decision still needs this.
"""

from ._common import UnsupportedStructure, adjacency, ring_cycle
from ._aromatic import (
    _atom_ring_membership,
    _edge_index,
    _ring_adjacency,
    _ring_path_order,
    find_aromatic_fused_core,
)

_TURN_BY_DIFF = {3: 0, 4: 1, 2: -1}


def _ring_turn_units(graph, atom_rings, ring_order, fusion_bonds_by_pair):
    """Signed 60-degree turn at each internal ring, in path order (one
    entry per ring strictly between the two chain endpoints)."""
    turns = []
    for pos in range(1, len(ring_order) - 1):
        ring_idx = ring_order[pos]
        cycle = ring_cycle(graph, list(atom_rings[ring_idx]))
        left_edge = _edge_index(cycle, *fusion_bonds_by_pair[frozenset((ring_order[pos - 1], ring_idx))])
        right_edge = _edge_index(cycle, *fusion_bonds_by_pair[frozenset((ring_idx, ring_order[pos + 1]))])
        diff = (left_edge - right_edge) % 6
        if diff not in _TURN_BY_DIFF:
            raise UnsupportedStructure(
                "an ortho-fused ring junction with an unexpected fusion-bond "
                "geometry is not supported (see P-25.3.2.3.3)"
            )
        turns.append(_TURN_BY_DIFF[diff])
    return turns


def max_rings_in_horizontal_row(turns, n):
    """P-25.3.2.3.3 criterion (a) alone, for a catacondensed non-branching
    chain of `n` rings. `turns`: this chain's `_ring_turn_units` result
    (length n - 2). Returns the largest number of rings that lie on one of
    the 3 possible horizontal-row axes."""
    if n <= 2:
        return n
    bond_axis = [0] * (n - 1)
    for k in range(1, n - 1):
        bond_axis[k] = bond_axis[k - 1] + turns[k - 1]
    best = 0
    for residue in range(3):
        count = 0
        for i in range(n):
            in_bond = bond_axis[i - 1] if i - 1 >= 0 else None
            out_bond = bond_axis[i] if i <= n - 2 else None
            if (in_bond is not None and in_bond % 3 == residue) or (
                out_bond is not None and out_bond % 3 == residue
            ):
                count += 1
        best = max(best, count)
    return best


def count_rings_in_horizontal_row(mol) -> int:
    """P-25.3.2.3.3 criterion (a) for a plain all-carbon aromatic mancude
    ring system: the number of rings in the orientation's horizontal row.
    Raises `UnsupportedStructure` for anything other than a catacondensed,
    non-branching, non-peri-fused chain (see module docstring)."""
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
    ring_order = _ring_path_order(adj, n)
    graph = adjacency(mol)
    turns = _ring_turn_units(graph, atom_rings, ring_order, fusion_bonds_by_pair)
    return max_rings_in_horizontal_row(turns, n)
