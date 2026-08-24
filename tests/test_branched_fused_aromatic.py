import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_triphenylene():
    # C18H12, ring-adjacency degree sequence [1, 1, 1, 3] confirmed via
    # RDKit (one ring fused to the other three, which are not fused to
    # each other) -- the shape _aromatic.py's own _ring_path_order names
    # as its motivating branched-topology example.
    assert smiles_to_iupac("c1ccc2c(c1)c1ccccc1c1ccccc21") == "triphenylene"


def test_substituted_triphenylene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)c1ccccc1c1ccccc21")


def test_unrelated_branched_shape_raises():
    # a different branched-fusion hydrocarbon (verified via RDKit: a
    # degree-3 ring-adjacency node, but not triphenylene's exact skeleton)
    # must still fall through to the ordinary "not supported" path, not
    # accidentally match triphenylene.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc2c(c1)c1ccc3ccccc3c1c1ccccc21")
