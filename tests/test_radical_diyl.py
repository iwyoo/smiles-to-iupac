import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-71.2.3 worked examples: both radical carbons are
        # chain termini.
        ("[CH2][CH2]", "ethane-1,2-diyl"),
        ("[CH2]C[CH2]", "propane-1,3-diyl"),
        ("[CH2]CC[CH2]", "butane-1,4-diyl"),
        ("[CH2]CCC[CH2]", "pentane-1,5-diyl"),
        # An interior (non-terminus) radical position combined with a
        # terminus one.
        ("[CH2][CH]C", "propane-1,2-diyl"),
        ("[CH2]C[CH]C", "butane-1,3-diyl"),
        # Two interior positions, neither a chain terminus.
        ("C[CH][CH]C", "butane-2,3-diyl"),
    ],
)
def test_chain_diradical_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[CH]1[CH]CCCC1", "cyclohexane-1,2-diyl"),
        ("[CH]1C[CH]CCC1", "cyclohexane-1,3-diyl"),
        ("[CH]1CC[CH]CC1", "cyclohexane-1,4-diyl"),
        ("[CH]1[CH]CCC1", "cyclopentane-1,2-diyl"),
    ],
)
def test_ring_diradical_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diradical_branch_point_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C([CH2])C")


def test_diradical_mixed_valence_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C[CH]")


def test_three_radical_centers_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2][CH][CH2]")


def test_diradical_substituted_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[CH]1[CH]CCC1")
