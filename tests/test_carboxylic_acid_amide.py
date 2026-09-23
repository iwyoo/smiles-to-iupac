import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_worked_example_10_amino_10_oxodecanoic_acid():
    # Blue Book P-66.1.1.3.3's own worked example.
    assert smiles_to_iupac("OC(=O)CCCCCCCCC(=O)N") == "10-amino-10-oxodecanoic acid"


def test_worked_example_malonamic_acid():
    # malonamic acid, PubChem CID 75367: '3-amino-3-oxopropanoic acid'.
    assert smiles_to_iupac("OC(=O)CC(=O)N") == "3-amino-3-oxopropanoic acid"


def test_amide_locant_from_other_end():
    # PubChem CID 12522: '4-amino-4-oxobutanoic acid'.
    assert smiles_to_iupac("OC(=O)CCC(=O)N") == "4-amino-4-oxobutanoic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("OC(=O)C(Cl)C(=O)N") == "3-amino-2-chloro-3-oxopropanoic acid"


def test_plain_carboxylic_acid_still_works():
    assert smiles_to_iupac("CCC(=O)O") == "propanoic acid"


def test_plain_amide_still_works():
    assert smiles_to_iupac("CCC(=O)N") == "propanamide"


def test_n_substituted_amide_raises():
    # An N-alkyl-substituted amide alongside the acid is out of scope for
    # this narrow slice (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CC(=O)NC")


def test_amide_on_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCCCC1C(=O)N")
