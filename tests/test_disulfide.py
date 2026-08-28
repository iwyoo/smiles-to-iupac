import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_disulfide():
    # PubChem PUG REST CID 12232, auto-generated name matches exactly --
    # note the enclosing parens (unlike the plain 'methylsulfanylmethane'
    # of the corresponding sulfide), same as `_diselenide.py`.
    assert smiles_to_iupac("CSSC") == "(methyldisulfanyl)methane"


def test_diethyl_disulfide():
    # PubChem PUG REST CID 8077, auto-generated name matches exactly.
    assert smiles_to_iupac("CCSSCC") == "(ethyldisulfanyl)ethane"


def test_methyl_ethyl_disulfide():
    # A 2-carbon parent omits the locant even though the sole substituent
    # is compound (parenthesized). PubChem PUG REST CID 123388,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CSSCC") == "(methyldisulfanyl)ethane"


def test_methyl_propyl_disulfide():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 16592,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCCSSC") == "1-(methyldisulfanyl)propane"


def test_branched_disulfanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)SSC(C)C")


def test_trisulfur_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CSSSC")
