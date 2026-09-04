import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_selenocyanate():
    # PubChem structure match: "methyl selenocyanate".
    assert smiles_to_iupac("C[Se]C#N") == "methyl selenocyanate"


def test_ethyl_selenocyanate():
    # PubChem structure match: "ethyl selenocyanate".
    assert smiles_to_iupac("CC[Se]C#N") == "ethyl selenocyanate"


def test_propyl_selenocyanate():
    # PubChem structure match: "propyl selenocyanate".
    assert smiles_to_iupac("CCC[Se]C#N") == "propyl selenocyanate"


def test_branched_r():
    # PubChem-verified: CID 13496974.
    assert smiles_to_iupac("CC(C)[Se]C#N") == "propan-2-yl selenocyanate"


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Se]C#N")


def test_cyclic_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1[Se]C#N")


def test_two_selenocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#C[Se]CC[Se]C#N")


def test_thiocyanate_not_confused_with_selenocyanate():
    assert smiles_to_iupac("CSC#N") == "methyl thiocyanate"


def test_cyanate_not_confused_with_selenocyanate():
    assert smiles_to_iupac("COC#N") == "methyl cyanate"


def test_phenyl_selenocyanate():
    # A plain, unsubstituted benzene ring bonded directly to the
    # selenocyanate selenium, cross-checked against PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1[Se]C#N") == "phenyl selenocyanate"  # CID 555340


def test_phenyl_selenocyanate_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C[Se]C#N")
