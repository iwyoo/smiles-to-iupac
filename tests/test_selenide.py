import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_selenide():
    # PubChem PUG REST CID 11648, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Se]C") == "methylselanylmethane"


def test_diethyl_selenide():
    # PubChem PUG REST CID 61173, auto-generated name matches exactly.
    assert smiles_to_iupac("CC[Se]CC") == "ethylselanylethane"


def test_ethyl_methyl_selenide():
    # The shorter (methyl) side becomes the substituent, the longer
    # (ethyl) side the parent. PubChem PUG REST CID 12248622,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("C[Se]CC") == "methylselanylethane"


def test_methyl_propyl_selenide():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 15932888,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[Se]C") == "1-methylselanylpropane"


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'selanyl' outside. Structure verified against PubChem: CID 11321295
    # ("1-propan-2-ylselanylbutane" -- PubChem elides the parentheses this
    # project's usual compound-substituent formatting keeps).
    assert smiles_to_iupac("CCCC[Se]C(C)C") == "1-(propan-2-yl)selanylbutane"


def test_both_sides_branched_and_tied_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se]C(C)C")
