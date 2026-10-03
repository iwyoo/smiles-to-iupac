import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Skeletal-carbon isotope (including the carboxyl carbon itself):
        # placed at the front of the whole name, per the Blue Book's own
        # '(1-14C)pentan(3H)oic acid (PIN)' worked example (P8.txt line 327).
        ("[13CH3]C(=O)O", "(2-13C)ethanoic acid"),
        ("[13CH3]CC(=O)O", "(3-13C)propanoic acid"),
        # Carboxyl-oxygen isotope: placed immediately before "oic acid",
        # after the suffix locant, matching the identical carbonyl-isotope
        # placement `_isotope_ketone.py` uses.
        ("CC(=[18O])O", "ethan(18O)oic acid"),
        ("CC(=O)[18OH]", "ethan(18O)oic acid"),
        # Both together on one name: two separate parenthetical groups
        # (front + before-suffix), exactly the shape the pentanoic-acid
        # worked example confirms.
        ("[14CH3]CCC(=[18O])O", "(4-14C)butan(18O)oic acid"),
    ],
)
def test_isotope_carboxylic_acid_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unlabeled_carboxylic_acid_unaffected():
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_isotope_ketone_unaffected_by_carboxylic_acid_routing():
    assert smiles_to_iupac("CC(=[18O])C") == "propan-2-(18O)one"


def test_ring_alongside_isotope_carboxylic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH2](C(=O)O)C1CCCCC1")


def test_chain_unsaturation_alongside_isotope_carboxylic_acid():
    assert smiles_to_iupac("[13CH3]C(=O)OC=C") == "ethenyl (2-13C)ethanoate"


def test_specified_stereocenter_alongside_isotope_carboxylic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3][C@H](Cl)C(=O)O")


def test_mixed_carbon_isotope_nuclides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3][14CH2]C(=O)O")
