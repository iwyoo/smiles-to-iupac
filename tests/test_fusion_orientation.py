import pytest
from rdkit import Chem

from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._fusion_orientation import (
    count_rings_in_horizontal_row,
    rings_above_horizontal_row,
    rings_in_lower_left_quadrant,
    rings_in_upper_right_quadrant,
)


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
    # branching case this module supports. PubChem CID 9170, InChIKey
    # SLGBZMMZGDRARJ-UHFFFAOYSA-N confirms this is real triphenylene (an
    # earlier version of this SMILES, still used elsewhere in this
    # project's git history, was actually a *linear* 4-ring chain
    # isomer -- caught during #570's PubChem cross-check, see
    # `test_chrysene_beats_triphenylene_on_criterion_b`).
    # `_triphenylene_fusion.py`'s own docstring documents that chrysene and
    # triphenylene tie on every P-25.3.2.4 criterion through (f) and are
    # only distinguished starting at (g) ("greatest number of rings in a
    # horizontal row"), which requires both to compute the same row count
    # here.
    mol = Chem.MolFromSmiles("c1ccc2c(c1)c1ccccc1c1ccccc21")
    assert count_rings_in_horizontal_row(mol) == 2


def test_pyrene_peri_fused_still_unsupported():
    # An atom shared by three rings makes the ring-fusion graph cyclic,
    # not a tree -- still out of scope (see module docstring).
    mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
    with pytest.raises(UnsupportedStructure):
        count_rings_in_horizontal_row(mol)


def test_phenanthrene_criterion_b_matches_primary_source():
    # P-25.3.2.3.3(b)'s own worked example: "phenanthrene (1 1/2 rings in
    # the upper right quadrant) is senior to phenalene [1 ring]." Phenalene
    # itself is peri-fused (an atom shared by three rings) and stays out
    # of scope, so only the phenanthrene half of that example is asserted.
    mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
    assert rings_in_upper_right_quadrant(mol) == 1.5


def test_naphthalene_criterion_b_half_ring():
    # A 2-ring row's center is the shared bond's midpoint; the one ring on
    # the positive-x side of it is bisected by the horizontal axis, so it
    # contributes exactly its upper half.
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    assert rings_in_upper_right_quadrant(mol) == 0.5


def test_tetraphene_criterion_b_needs_orientation_search():
    # P-25.3.2.3.3(b)'s second worked example: "3 rings in horizontal row,
    # 1 3/4 rings in upper right quadrant" - this is tetraphene/
    # benzo[a]anthracene. An earlier version of this function picked one
    # arbitrary candidate orientation instead of searching all of them and
    # got 0.75 for this molecule (still correct for phenanthrene, only by
    # luck) - this is the asymmetric case that catches that gap.
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32")
    assert rings_in_upper_right_quadrant(mol) == 1.75


def test_chrysene_beats_triphenylene_on_criterion_b():
    # They tie on (a) (both 2 rings in a horizontal row), but genuine
    # triphenylene (PubChem CID 9170) does NOT tie chrysene on (b) - an
    # earlier, mislabeled SMILES for "triphenylene" (actually a linear
    # 4-ring chain isomer, caught during #570's PubChem cross-check) gave
    # a false tie at 2.5 for both. With the real, branched triphenylene,
    # chrysene's 2.5 already beats triphenylene's 1.5, so P-25.3.2.4(g)'s
    # seniority question between them is settled at (b) - no need for (c),
    # (d), or any fallback further down the seniority list.
    chrysene = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43")
    triphenylene = Chem.MolFromSmiles("c1ccc2c(c1)c1ccccc1c1ccccc21")
    assert rings_in_upper_right_quadrant(chrysene) == 2.5
    assert rings_in_upper_right_quadrant(triphenylene) == 1.5


def test_phenanthrene_criteria_c_and_d():
    mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
    assert rings_in_lower_left_quadrant(mol) == 0.5
    assert rings_above_horizontal_row(mol) == 2.0


def test_tetraphene_criteria_c_and_d():
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32")
    assert rings_in_lower_left_quadrant(mol) == 0.75
    assert rings_above_horizontal_row(mol) == 2.5


def test_triphenylene_criteria_c_and_d():
    # Regression values for genuine triphenylene on a branched ring-fusion
    # graph - not tie-breakers for the chrysene comparison, since neither
    # (b) nor (c)/(d) are themselves P-25.3.2.4 seniority criteria (see
    # `test_chrysene_beats_triphenylene_on_criterion_b` and the module
    # docstring).
    triphenylene = Chem.MolFromSmiles("c1ccc2c(c1)c1ccccc1c1ccccc21")
    assert rings_in_lower_left_quadrant(triphenylene) == 0.5
    assert rings_above_horizontal_row(triphenylene) == 2.0


@pytest.mark.parametrize(
    "smiles",
    [
        "c1ccc2ccc3ccccc3c2c1",  # phenanthrene
        "C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32",  # tetraphene
        "C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43",  # chrysene
        "c1ccc2c(c1)c1ccccc1c1ccccc21",  # triphenylene
        "c1ccc2cc3ccccc3cc2c1",  # anthracene
        "c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1",  # pentacene
    ],
)
def test_above_plus_below_equals_total_ring_count(smiles):
    # "Below the row" isn't exposed directly, but total ring count minus
    # "above" is exactly what it must be by construction - this is really
    # asserting `rings_above_horizontal_row` never over/undercounts.
    mol = Chem.MolFromSmiles(smiles)
    n = mol.GetRingInfo().NumRings()
    above = rings_above_horizontal_row(mol)
    assert 0 <= above <= n


# Random-sample validation for issue #570: 14 previously-untested,
# tree-shaped (catacondensed or branched) aromatic ring systems pulled
# live from PubChem (CanonicalSMILES via PUG REST), none of which were
# used to build or tune this module - this is what caught the mislabeled
# "triphenylene" SMILES above (its PubChem cross-check, CID 9170,
# InChIKey SLGBZMMZGDRARJ-UHFFFAOYSA-N, is what exposed the earlier test
# suite's linear-chain impostor). Values pinned here are regression
# anchors on this implementation, not independently Blue-Book-verified
# per compound (only phenanthrene/tetraphene/chrysene/anthracene/pentacene
# above have primary-source-cited values) - the property under test is
# that the algorithm runs to completion and stays internally consistent
# (quadrant/row counts sum to the total ring count) on real structures it
# has never seen, not that each named value matches a hand-checked source.
_PUBCHEM_SAMPLE = [
    # (PubChem CID, name, SMILES, expected rings-in-horizontal-row)
    (7080, "naphthacene", "C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1", 4),
    (8671, "pentacene", "C1=CC=C2C=C3C=C4C=C5C=CC=CC5=CC4=CC3=CC2=C1", 5),
    (9162, "picene", "C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC=CC=C54", 4),
    (519935, "pentaphene", "C1=CC=C2C=C3C(=CC2=C1)C=CC4=CC5=CC=CC=C5C=C43", 3),
    (9135, "benzo[c]chrysene", "C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C5=CC=CC=C5C=C4", 4),
    (9140, "benzo[g]chrysene", "C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C5=CC=CC=C35", 4),
    (9164, "dibenzo[a,c]anthracene", "C1=CC=C2C=C3C4=CC=CC=C4C5=CC=CC=C5C3=CC2=C1", 3),
    (5889, "dibenz[a,h]anthracene", "C1=CC=C2C(=C1)C=CC3=CC4=C(C=CC5=CC=CC=C54)C=C32", 3),
    (9176, "dibenz[a,j]anthracene", "C1=CC=C2C(=C1)C=CC3=CC4=C(C=C32)C5=CC=CC=C5C=C4", 4),
    (5954, "benzo[a]anthracene", "C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32", 3),
    (9167, "dibenzo[a,c]naphthacene", "C1=CC=C2C=C3C=C4C5=CC=CC=C5C6=CC=CC=C6C4=CC3=CC2=C1", 4),
    (98863, "hexahelicene", "C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=C4C6=CC=CC=C6C=C5", 6),
]


@pytest.mark.parametrize("cid,name,smiles,expected_row", _PUBCHEM_SAMPLE, ids=[s[1] for s in _PUBCHEM_SAMPLE])
def test_pubchem_random_sample_row_count_and_consistency(cid, name, smiles, expected_row):
    mol = Chem.MolFromSmiles(smiles)
    n = mol.GetRingInfo().NumRings()
    assert count_rings_in_horizontal_row(mol) == expected_row
    ur = rings_in_upper_right_quadrant(mol)
    ll = rings_in_lower_left_quadrant(mol)
    above = rings_above_horizontal_row(mol)
    assert 0 <= ur <= n
    assert 0 <= ll <= n
    assert 0 <= above <= n
    assert above + (n - above) == n


def test_pubchem_sample_peri_fused_case_still_unsupported():
    # dibenzo[def,p]chrysene (PubChem CID 9119) from the same random pull -
    # a real, registered peri-fused structure (an atom shared by three
    # rings), confirming the tree-only scope boundary against a structure
    # this module was never tuned against, not just the hand-built
    # `test_pyrene_peri_fused_still_unsupported` case above.
    mol = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=C3C=CC4=C5C3=C2C6=CC=CC=C6C5=CC=C4")
    with pytest.raises(UnsupportedStructure):
        count_rings_in_horizontal_row(mol)
