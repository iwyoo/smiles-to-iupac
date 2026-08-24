import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_worked_example_malonaldehydic_acid():
    # malonaldehydic acid / malonic semialdehyde (PubChem CID 868): IUPAC
    # name '3-oxopropanoic acid', confirmed via PubChem.
    assert smiles_to_iupac("O=CCC(=O)O") == "3-oxopropanoic acid"


def test_aldehyde_locant_from_other_end():
    assert smiles_to_iupac("OC(=O)CCC=O") == "4-oxobutanoic acid"


def test_branch_substituent():
    assert smiles_to_iupac("O=CC(C)C(=O)O") == "2-methyl-3-oxopropanoic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("O=CC(Cl)C(=O)O") == "2-chloro-3-oxopropanoic acid"


def test_plain_carboxylic_acid_still_works():
    assert smiles_to_iupac("CCC(=O)O") == "propanoic acid"


def test_plain_aldehyde_still_works():
    assert smiles_to_iupac("CCC=O") == "propanal"


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCC(=O)CC(=O)O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC1=O")


def test_hydroxyl_coexistence_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=O)CC(=O)O")
