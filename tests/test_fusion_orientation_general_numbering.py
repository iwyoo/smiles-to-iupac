from fractions import Fraction

from rdkit import Chem

from smiles_to_iupac._fusion_orientation import _prepare
from smiles_to_iupac._fusion_orientation_general import starting_ring_general


def _general_direction(mol):
    _atom_rings, adj, _fusion_bonds_by_pair, direction, n = _prepare(mol)
    return adj, {pair: Fraction(d, 6) % 1 for pair, d in direction.items()}, n


def test_anthracene_starting_ring_is_a_terminal_ring():
    # Straight 3-ring chain: all three ring centers lie on the same row,
    # so P-25.3.3.1.1's "uppermost" tie-break falls through to "farthest
    # right" -- one specific terminal ring, never the symmetric middle one.
    mol = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
    adj, direction, n = _general_direction(mol)
    winners = starting_ring_general(adj, direction, n)
    assert len(winners) == 1
    assert len(adj[winners[0]]) == 1


def test_phenanthrene_starting_ring_is_unique():
    # Angular 3-ring chain: geometry alone (no heteroatom tie-break
    # needed) picks exactly one ring as uppermost-then-rightmost.
    mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
    adj, direction, n = _general_direction(mol)
    winners = starting_ring_general(adj, direction, n)
    assert len(winners) == 1


def test_pentagon_between_two_hexagons_ties_like_anthracene():
    # Same straight-chain topology as anthracene, but the middle ring is a
    # pentagon (5) instead of a hexagon -- P-25.3.2.3.1's own permitted
    # shape for this size still keeps the chain's mirror symmetry, so the
    # starting-ring tie-break behaves the same way as the all-hexagon case.
    ring_sizes = {0: 6, 1: 5, 2: 6}
    adj = {0: {1}, 1: {0, 2}, 2: {1}}
    edge_index_of = {(0, 1): 0, (1, 0): 0, (1, 2): 2, (2, 1): 0}
    from smiles_to_iupac._fusion_orientation_general import assign_bond_directions_general

    direction = assign_bond_directions_general(adj, edge_index_of, ring_sizes, 3)
    winners = starting_ring_general(adj, direction, 3)
    assert len(winners) == 1
    assert len(adj[winners[0]]) == 1
