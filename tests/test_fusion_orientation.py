import pytest
from rdkit import Chem

from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._fusion_orientation import count_rings_in_horizontal_row


def test_anthracene_all_three_in_row():
    # P-25.3.2.3.3's own example: "anthracene is senior to phenanthrene"
    # (3 rings in a row vs 2).
    mol = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
    assert count_rings_in_horizontal_row(mol) == 3


def test_phenanthrene_two_in_row():
    mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
    assert count_rings_in_horizontal_row(mol) == 2


def test_tetraphene_three_in_row():
    # P-25.3.2.4's own worked example: "tetraphene (3 rings in horizontal
    # row) preferred to chrysene or pyrene (2 rings in horizontal row)".
    # Tetraphene is benzo[a]anthracene (PubChem CID 5954).
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32")
    assert count_rings_in_horizontal_row(mol) == 3


def test_chrysene_two_in_row():
    # PubChem CID 9171.
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43")
    assert count_rings_in_horizontal_row(mol) == 2


def test_benzo_c_phenanthrene_distinguished_from_chrysene():
    # PubChem CID 9136 -- collapses to the same ('bent', 'bent') shape as
    # chrysene under `_aromatic.py`'s `_classify_shape`, but its two bends
    # turn in opposite directions (diffs [4, 2]) rather than the same
    # direction twice (chrysene's [4, 4]). This is exactly why this module
    # recomputes signed turns instead of consuming `_classify_shape`'s
    # already-collapsed labels -- no primary-source value is asserted here
    # (only the Blue Book's own four worked examples are used as ground
    # truth elsewhere in this file), but a correct model must not treat
    # these two structurally distinct compounds identically by accident.
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C=C3")
    chrysene = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43")
    assert count_rings_in_horizontal_row(mol) != count_rings_in_horizontal_row(chrysene)


def test_naphthalene_trivially_two():
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    assert count_rings_in_horizontal_row(mol) == 2


def test_benzene_trivially_one():
    mol = Chem.MolFromSmiles("c1ccccc1")
    assert count_rings_in_horizontal_row(mol) == 1


def test_straight_tetracene_all_four_in_row():
    mol = Chem.MolFromSmiles("c1ccc2cc3cc4ccccc4cc3cc2c1")
    assert count_rings_in_horizontal_row(mol) == 4


def test_straight_pentacene_all_five_in_row():
    mol = Chem.MolFromSmiles("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1")
    assert count_rings_in_horizontal_row(mol) == 5


def test_triphenylene_branching_two_in_row():
    # Triphenylene's ring-fusion graph is branched (one central ring
    # ortho-fused to three others), not a simple chain -- the first
    # branching case this module supports. `_triphenylene_fusion.py`'s own
    # docstring documents that chrysene and triphenylene tie on every
    # criterion through (f) and are only distinguished starting at (g)
    # ("greatest number of rings in a horizontal row"), which requires
    # both to compute the same row count here.
    mol = Chem.MolFromSmiles("c1ccc2c(c1)ccc1c2ccc2ccccc21")
    assert count_rings_in_horizontal_row(mol) == 2


def test_pyrene_peri_fused_still_unsupported():
    # An atom shared by three rings makes the ring-fusion graph cyclic,
    # not a tree -- still out of scope (see module docstring).
    mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
    with pytest.raises(UnsupportedStructure):
        count_rings_in_horizontal_row(mol)
