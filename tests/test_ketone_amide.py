import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_worked_example_acetoacetamide():
    # acetoacetamide (PubChem CID 80077, CAS 5977-14-0): IUPAC name
    # '3-oxobutanamide', confirmed via PubChem.
    assert smiles_to_iupac("CC(=O)CC(N)=O") == "3-oxobutanamide"


def test_ketone_locant_from_other_end():
    assert smiles_to_iupac("NC(=O)CCC(=O)C") == "4-oxopentanamide"


def test_two_ketones():
    assert smiles_to_iupac("CC(=O)CC(=O)CC(N)=O") == "3,5-dioxohexanamide"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(=O)C(Cl)C(N)=O") == "2-chloro-3-oxobutanamide"


def test_plain_amide_still_works():
    assert smiles_to_iupac("CCC(N)=O") == "propanamide"


def test_plain_ketone_still_works():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)CC(N)=O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)C1CCC1=O")


def test_hydroxyl_coexistence_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=O)CC(N)=O")


def test_n_substituted_amide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(=O)CC(=O)C")
