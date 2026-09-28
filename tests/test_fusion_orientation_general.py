import math
from fractions import Fraction

from rdkit import Chem

from smiles_to_iupac._aromatic import _edge_index, _ring_adjacency, find_aromatic_fused_core
from smiles_to_iupac._common import UnsupportedStructure, adjacency, ring_cycle
from smiles_to_iupac._fusion_orientation import (
    count_rings_in_horizontal_row,
    rings_above_horizontal_row,
    rings_in_lower_left_quadrant,
    rings_in_upper_right_quadrant,
)
from smiles_to_iupac._fusion_orientation_general import (
    assign_bond_directions_general,
    best_orientation_general,
    max_rings_in_horizontal_row_general,
)


def _general_inputs_from_mol(mol):
    """Build (adj, edge_index_of, ring_sizes, n) for the general algorithm
    from a real hexagon-only molecule, using the same `_aromatic.py`
    extraction the production hexagon-only module uses -- lets a real
    PubChem structure directly cross-check the general algorithm against
    `_fusion_orientation.py`'s own hexagon-only production results."""
    atom_rings, ring_atom_sets, fusion_bond_idxs = find_aromatic_fused_core(mol)
    n = len(atom_rings)
    adj, fusion_bonds_by_pair = _ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol)
    ring_sizes = {i: len(atom_rings[i]) for i in range(n)}
    graph = adjacency(mol)
    edge_index_of = {}
    for pair, (a, b) in fusion_bonds_by_pair.items():
        i, j = tuple(pair)
        cycle_i = ring_cycle(graph, list(atom_rings[i]))
        cycle_j = ring_cycle(graph, list(atom_rings[j]))
        edge_index_of[(i, j)] = _edge_index(cycle_i, a, b)
        edge_index_of[(j, i)] = _edge_index(cycle_j, a, b)
    return adj, edge_index_of, ring_sizes, n


def _general_results(mol):
    adj, edge_index_of, ring_sizes, n = _general_inputs_from_mol(mol)
    direction = assign_bond_directions_general(adj, edge_index_of, ring_sizes, n)
    row = max_rings_in_horizontal_row_general(adj, direction, n)
    _residue, _reflect, upper_right, lower_left, above = best_orientation_general(adj, direction, n)
    return row, upper_right, lower_left, above


def test_hexagon_only_reduction_matches_existing_module_anthracene():
    # PubChem CID 8418. The general (Fraction-based, arbitrary-ring-size)
    # algorithm must reproduce the production hexagon-only module's
    # result exactly when every ring is N=6.
    mol = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
    row, _ur, _ll, _above = _general_results(mol)
    assert row == count_rings_in_horizontal_row(mol) == 3


def test_hexagon_only_reduction_matches_existing_module_phenanthrene():
    mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
    row, ur, _ll, _above = _general_results(mol)
    assert row == count_rings_in_horizontal_row(mol) == 2
    assert ur == rings_in_upper_right_quadrant(mol) == 1.5


def test_hexagon_only_reduction_matches_existing_module_tetraphene():
    # PubChem CID 5954 (benzo[a]anthracene).
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32")
    row, ur, _ll, _above = _general_results(mol)
    assert row == count_rings_in_horizontal_row(mol) == 3
    assert ur == rings_in_upper_right_quadrant(mol) == 1.75


def test_hexagon_only_reduction_matches_existing_module_chrysene():
    # PubChem CID 9171.
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43")
    row, ur, ll, above = _general_results(mol)
    assert row == count_rings_in_horizontal_row(mol) == 2
    assert ur == rings_in_upper_right_quadrant(mol)
    assert ll == rings_in_lower_left_quadrant(mol)
    assert above == rings_above_horizontal_row(mol)


def _synthetic_chain(sizes, edge_pair_at_middle):
    """A 3-ring linear chain ring0-ring1-ring2 where ring1 (the middle
    ring, of size `sizes[1]`) connects to ring0 via edge index
    `edge_pair_at_middle[0]` and to ring2 via `edge_pair_at_middle[1]` (in
    ring1's own cyclic bond order) -- ring0/ring2's own edge indices don't
    affect any result here (a leaf ring's `ref_edge` is never used to
    compute anything further)."""
    adj = {0: {1}, 1: {0, 2}, 2: {1}}
    ring_sizes = {0: sizes[0], 1: sizes[1], 2: sizes[2]}
    edge_index_of = {
        (0, 1): 0,
        (1, 0): edge_pair_at_middle[0],
        (1, 2): edge_pair_at_middle[1],
        (2, 1): 0,
    }
    return adj, edge_index_of, ring_sizes, 3


def test_hexagon_chain_straight_matches_anthracene_row_three():
    # Hand derivation: with ring0's fusion bond at edge 0 and ring1's own
    # opposite edge at k=3 (ceil(6/2)), diff=(0-3)%6=3, centered=3-3=0 --
    # direction[{1,2}] == direction[{0,1}] == 0, same axis, row=3,
    # matching anthracene's real row=3 (P-25.3.2.3.3's own "anthracene is
    # senior to phenanthrene" example).
    adj, edge_index_of, ring_sizes, n = _synthetic_chain((6, 6, 6), (0, 3))
    direction = assign_bond_directions_general(adj, edge_index_of, ring_sizes, n)
    assert direction[frozenset((1, 2))] == Fraction(0)
    assert max_rings_in_horizontal_row_general(adj, direction, n) == 3


def test_hexagon_chain_bent_matches_phenanthrene_row_two():
    # Hand derivation: k=2 (not the unique straight-through k=3): diff=
    # (0-2)%6=4, centered=4-3=1, direction[{1,2}]=Fraction(1,6) -- a
    # DIFFERENT axis than direction[{0,1}]=0 -- so only 2 rings share any
    # one axis, matching phenanthrene's real row=2.
    adj, edge_index_of, ring_sizes, n = _synthetic_chain((6, 6, 6), (0, 2))
    direction = assign_bond_directions_general(adj, edge_index_of, ring_sizes, n)
    assert direction[frozenset((1, 2))] == Fraction(1, 6)
    assert max_rings_in_horizontal_row_general(adj, direction, n) == 2


def test_exactly_one_edge_gap_per_ring_size_gives_straight_row_three():
    # For every permitted P-25.3.2.3.1 ring size N=3..8, there is exactly
    # one edge-index gap k (between the middle ring's two fusion bonds)
    # that achieves criterion (a)'s row=3 "straight-through" alignment --
    # always k = ceil(N/2), the ring's own closest analogue to a
    # hexagon's antipodal edge pair (verified here for every k, not
    # assumed): N=6 reproduces the real anthracene/phenanthrene distinction
    # above; N=3,4,5,7,8 have no real-molecule ground truth available (see
    # module docstring), so these are self-consistency checks of the same
    # formula, not independently verified against an external source.
    for N in (3, 4, 5, 6, 7, 8):
        expected_straight_k = math.ceil(N / 2)
        straight_ks = []
        for k in range(1, N):
            adj, edge_index_of, ring_sizes, n = _synthetic_chain((6, N, 6), (0, k))
            direction = assign_bond_directions_general(adj, edge_index_of, ring_sizes, n)
            if max_rings_in_horizontal_row_general(adj, direction, n) == 3:
                straight_ks.append(k)
        assert straight_ks == [expected_straight_k], N


def test_disconnected_ring_fusion_graph_rejected():
    adj = {0: set(), 1: set()}
    edge_index_of = {}
    ring_sizes = {0: 6, 1: 6}
    try:
        assign_bond_directions_general(adj, edge_index_of, ring_sizes, 2)
        assert False, "expected UnsupportedStructure"
    except UnsupportedStructure:
        pass
