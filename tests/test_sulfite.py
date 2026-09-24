import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-67.1.3.2), symmetric dialkyl/diaryl
        # sulfite esters -- verified against each structure's own
        # PubChem-registered IUPAC name.
        ("COS(=O)OC", "dimethyl sulfite"),  # CID 69223
        ("CCOS(=O)OCC", "diethyl sulfite"),  # CID 12197
        ("CCCOS(=O)OCCC", "dipropyl sulfite"),  # CID 136434
        ("c1ccc(OS(=O)Oc2ccccc2)cc1", "diphenyl sulfite"),  # CID 11390667
    ],
)
def test_sulfite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_alkyl_sulfite_ester_name():
    # Two distinct R groups: two separate un-prefixed words, alphabetic
    # order (same shape `_sulfate.py`'s mixed-alkyl case tests).
    assert smiles_to_iupac("CCOS(=O)OC") == "ethyl methyl sulfite"


def test_partial_hydrogen_sulfite_ester_name():
    # Real PubChem structure (P-67.1.3.2), one R group + one remaining
    # S-OH.
    assert smiles_to_iupac("COS(=O)O") == "methyl hydrogen sulfite"  # CID 358915


def test_sulfurous_acid_itself_raises():
    # Both S-O positions are plain hydroxyl (no R group at all) -- the
    # parent acid itself, not an ester, still out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)O")


def test_sulfate_ester_unaffected():
    # A second S=O bond routes to `_sulfate.py` instead, unchanged.
    assert smiles_to_iupac("COS(=O)(=O)OC") == "dimethyl sulfate"
