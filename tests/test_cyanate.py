import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_cyanate():
    # PubChem structure match: "methyl cyanate".
    assert smiles_to_iupac("COC#N") == "methyl cyanate"


def test_ethyl_cyanate():
    # PubChem structure match: "ethyl cyanate".
    assert smiles_to_iupac("CCOC#N") == "ethyl cyanate"


def test_propyl_cyanate():
    # PubChem structure match: "propyl cyanate".
    assert smiles_to_iupac("CCCOC#N") == "propyl cyanate"


def test_branched_r():
    # PubChem-verified: CID 550695.
    assert smiles_to_iupac("CC(C)OC#N") == "propan-2-yl cyanate"


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=COC#N")


def test_cyclic_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1OC#N")


def test_two_cyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#COCCOC#N")


def test_ether_not_confused_with_cyanate():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_thiocyanate_not_confused_with_cyanate():
    assert smiles_to_iupac("CSC#N") == "methyl thiocyanate"


def test_phenyl_cyanate():
    # A plain, unsubstituted benzene ring bonded directly to the cyanate
    # oxygen, cross-checked against PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1OC#N") == "phenyl cyanate"  # CID 70740


def test_phenyl_cyanate_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1COC#N")


def test_phenyl_cyanate_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1OC#N")
