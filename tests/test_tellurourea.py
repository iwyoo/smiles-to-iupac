import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_tellurourea():
    # PubChem structure match: NC(=[Te])N is CID 14205311. Blue Book
    # P-66.1.6.1.3.1's rule text names 'thio', 'seleno', and 'telluro'
    # together as the three chalcogen analogue prefixes, so 'tellurourea'
    # is a PIN by the same mechanism as 'thiourea'/'selenourea', even
    # without a tellurourea-specific worked example.
    assert smiles_to_iupac("NC(=[Te])N") == "tellurourea"


def test_n_methyltellurourea():
    assert smiles_to_iupac("CNC(=[Te])N") == "N-methyltellurourea"


def test_n_n_dimethyltellurourea_same_nitrogen():
    assert smiles_to_iupac("CN(C)C(=[Te])N") == "N,N-dimethyltellurourea"


def test_n_ethyl_n_methyltellurourea_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=[Te])N") == "N-ethyl-N-methyltellurourea"


def test_n_n_prime_dimethyltellurourea_different_nitrogens():
    assert smiles_to_iupac("CNC(=[Te])NC") == "N,N'-dimethyltellurourea"


def test_different_substituents_on_different_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCNC(=[Te])NC")


def test_branched_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)NC(=[Te])N")


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=[Te])N")


def test_selenourea_not_confused_with_tellurourea():
    assert smiles_to_iupac("NC(=[Se])N") == "selenourea"


def test_urea_not_confused_with_tellurourea():
    assert smiles_to_iupac("NC(=O)N") == "urea"
