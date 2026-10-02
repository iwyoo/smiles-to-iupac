from fractions import Fraction

from rdkit import Chem

from smiles_to_iupac._fusion_orientation import _prepare
from smiles_to_iupac._fusion_orientation_general import starting_ring_general


def _general_direction(mol):
    _atom_rings, adj, _fusion_bonds_by_pair, direction, n = _prepare(mol)
    return adj, {pair: Fraction(d, 6) % 1 for pair, d in direction.items()}, n


def test_anthracene_starting_ring_is_a_terminal_ring():
    # Straight 3-ring chain: both terminal rings tie as "uppermost,
    # farthest right" (anthracene's mirror symmetry swaps them, and
    # criteria (a)-(d) never compare ring shapes) -- never the middle one.
    mol = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
    adj, direction, n = _general_direction(mol)
    winners = starting_ring_general(adj, direction, n)
    assert all(len(adj[w]) == 1 for w in winners)


def test_phenanthrene_starting_ring_excludes_no_terminal_ring():
    # Angular 3-ring chain: phenanthrene also has a mirror symmetry
    # swapping its two terminal rings, so both legitimately tie through
    # criteria (a)-(d) -- same reasoning as the anthracene case above.
    mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
    adj, direction, n = _general_direction(mol)
    winners = starting_ring_general(adj, direction, n)
    assert all(len(adj[w]) == 1 for w in winners)
    assert len(winners) >= 1


def test_pentagon_between_two_hexagons_excludes_middle_ring():
    # Same straight-chain topology as anthracene, but the middle ring is a
    # pentagon (5) instead of a hexagon -- P-25.3.2.3.1's own permitted
    # shape for this size still keeps the chain's mirror symmetry, so the
    # starting-ring tie-break excludes the middle ring the same way the
    # all-hexagon case does.
    ring_sizes = {0: 6, 1: 5, 2: 6}
    adj = {0: {1}, 1: {0, 2}, 2: {1}}
    edge_index_of = {(0, 1): 0, (1, 0): 0, (1, 2): 2, (2, 1): 0}
    from smiles_to_iupac._fusion_orientation_general import assign_bond_directions_general

    direction = assign_bond_directions_general(adj, edge_index_of, ring_sizes, 3)
    winners = starting_ring_general(adj, direction, 3)
    assert all(len(adj[w]) == 1 for w in winners)
