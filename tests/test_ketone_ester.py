import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_worked_example_methyl_acetoacetate():
    # methyl acetoacetate (PubChem CID 7757, CAS 105-45-3): IUPAC name
    # 'methyl 3-oxobutanoate', confirmed via PubChem/NIST WebBook.
    assert smiles_to_iupac("CC(=O)CC(=O)OC") == "methyl 3-oxobutanoate"


def test_ethyl_acetoacetate():
    # A single-axis variant of the worked example above (methyl -> ethyl
    # ester), same demotion mechanism.
    assert smiles_to_iupac("CCOC(=O)CC(=O)C") == "ethyl 3-oxobutanoate"


def test_two_ketones():
    assert smiles_to_iupac("CC(=O)CC(=O)CC(=O)OC") == "methyl 3,5-dioxohexanoate"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(=O)C(Cl)C(=O)OC") == "methyl 2-chloro-3-oxobutanoate"


def test_plain_ester_still_works():
    assert smiles_to_iupac("CCC(=O)OC") == "methyl propanoate"


def test_plain_ketone_still_works():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)CC(=O)OC")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)C1CCC1=O")


def test_hydroxyl_coexistence_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=O)CC(=O)OC")
