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


def test_quaternary_phosphonium_not_supported():
    # A quaternary phosphonium (four carbons) has no neutral phosphane
    # counterpart to derive its name from -- structurally confirmed on
    # PubChem as "tetramethylphosphanium", but out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+](C)(C)C")


def test_ammonium_not_confused_with_phosphonium():
    assert smiles_to_iupac("C[NH3+]") == "methanaminium"


def test_phosphane_not_confused_with_phosphonium():
    assert smiles_to_iupac("CP") == "methylphosphane"
