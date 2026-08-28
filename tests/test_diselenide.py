import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_diselenide():
    # PubChem PUG REST CID 23496, auto-generated name matches exactly --
    # note the enclosing parens (unlike the plain 'methylselanylmethane'
    # of the corresponding selenide), see module docstring.
    assert smiles_to_iupac("C[Se][Se]C") == "(methyldiselanyl)methane"


def test_diethyl_diselenide():
    # PubChem PUG REST CID 69405, auto-generated name matches exactly.
    assert smiles_to_iupac("CC[Se][Se]CC") == "(ethyldiselanyl)ethane"


def test_methyl_propyl_diselenide():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 85591054,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[Se][Se]C") == "1-(methyldiselanyl)propane"


def test_methyl_ethyl_diselenide():
    # A 2-carbon parent omits the locant even though the sole substituent
    # is compound (parenthesized) -- see module docstring. PubChem PUG
    # REST CID 129678518, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Se][Se]CC") == "(methyldiselanyl)ethane"


def test_branched_diselanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se][Se]C(C)C")


def test_triselenium_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Se][Se][Se]C")
