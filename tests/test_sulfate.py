import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-67.1.3.2), symmetric dialkyl/diaryl
        # sulfate esters -- verified against each structure's own
        # PubChem-registered IUPAC name.
        ("COS(=O)(=O)OC", "dimethyl sulfate"),  # CID 6497
        ("CCOS(=O)(=O)OCC", "diethyl sulfate"),  # CID 6163
        ("CCCOS(=O)(=O)OCCC", "dipropyl sulfate"),  # CID 69015
        ("CC(C)OS(=O)(=O)OC(C)C", "dipropan-2-yl sulfate"),  # CID 18096
        ("c1ccc(OS(=O)(=O)Oc2ccccc2)cc1", "diphenyl sulfate"),  # CID 14228015
        # A branch name starting with its own locant digit needs 'bis'
        # plus enclosing marks (same rule `_phosphate.py`'s 'tris' case
        # tests, one fewer R group here).
        ("ClCCOS(=O)(=O)OCCCl", "bis(2-chloroethyl) sulfate"),  # CID 79428
    ],
)
def test_sulfate_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_alkyl_sulfate_ester_name():
    # Two distinct R groups: two separate un-prefixed words, alphabetic
    # order. PubChem CID 69945.
    assert smiles_to_iupac("CCOS(=O)(=O)OC") == "ethyl methyl sulfate"


def test_partial_hydrogen_sulfate_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COS(=O)(=O)O")


def test_sulfite_ester_still_unsupported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COS(=O)OC")
