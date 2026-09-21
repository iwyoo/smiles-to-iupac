import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Carbonyl-oxygen isotope: placed immediately before "one", after
        # the suffix locant, by direct analogy with the Blue Book's own
        # '1-(aminomethyl)cyclopentan-1-(18O)ol (PIN)' (the oxygen is the
        # suffix-defining atom in both cases).
        ("CC(=[18O])C", "propan-2-(18O)one"),
        ("CCC(=[18O])C", "butan-2-(18O)one"),
        # Skeletal-carbon isotope: placed at the front of the whole name,
        # per P-82.2.5 and the Blue Book's own '1-phenyl(1,2-13C2)ethan-
        # 1-one (PIN)' worked example.
        ("[13CH3]C(=O)C", "(1-13C)propan-2-one"),
        ("CC(=O)[14CH3]", "(3-14C)propan-2-one"),
        # A halogen substituent coexisting with a carbonyl-oxygen isotope.
        ("ClCC(=[18O])C", "1-chloropropan-2-(18O)one"),
        # A plain (non-isotope-labeled) ketone with unrelated unsaturation
        # elsewhere is unaffected by this module's routing.
    ],
)
def test_isotope_ketone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_both_isotope_kinds_together_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=[18O])[13CH3]")


def test_unlabeled_ketone_unaffected():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


def test_ring_alongside_isotope_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[13C]1CCCCC1")


def test_chain_unsaturation_alongside_isotope_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3]C(=O)C=C")


def test_specified_stereocenter_alongside_isotope_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3]C(=O)[C@H](Cl)C")


def test_mixed_carbon_isotope_nuclides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3]C(=O)[14CH3]")
