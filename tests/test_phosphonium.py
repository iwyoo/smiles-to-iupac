import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_phosphanium():
    # PubChem structure match: "phosphanium" (PH4+).
    assert smiles_to_iupac("[PH4+]") == "phosphanium"


def test_methylphosphanium():
    # PubChem structure match: "methylphosphanium".
    assert smiles_to_iupac("C[PH3+]") == "methylphosphanium"


def test_dimethylphosphanium():
    # PubChem structure match: "dimethylphosphanium".
    assert smiles_to_iupac("C[PH2+]C") == "dimethylphosphanium"


def test_trimethylphosphanium():
    # PubChem structure match: "trimethylphosphanium".
    assert smiles_to_iupac("C[PH+](C)C") == "trimethylphosphanium"


def test_ethylphosphanium():
    assert smiles_to_iupac("CC[PH3+]") == "ethylphosphanium"


def test_tetramethylphosphanium():
    # A quaternary phosphonium (four carbons) has no neutral phosphane
    # counterpart, so it's named directly from its substituents. PubChem
    # structure match: "tetramethylphosphanium".
    assert smiles_to_iupac("C[P+](C)(C)C") == "tetramethylphosphanium"


def test_ethyl_trimethyl_phosphanium():
    # PubChem structure match: "ethyl(trimethyl)phosphanium" -- the
    # P-16.5.1.3.1 parenthesization rule for a mononuclear parent with
    # mixed distinct substituent counts.
    assert smiles_to_iupac("CC[P+](C)(C)C") == "ethyl(trimethyl)phosphanium"


def test_branched_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+](C)(C)C(C)C")


def test_ring_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+]1(C)CCCC1")


def test_halogen_substituted_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+](C)(C)Cl")


def test_ammonium_not_confused_with_phosphonium():
    assert smiles_to_iupac("C[NH3+]") == "methanaminium"


def test_phosphane_not_confused_with_phosphonium():
    assert smiles_to_iupac("CP") == "methylphosphane"
