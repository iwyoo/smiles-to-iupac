import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_pyrene():
    # cross-checked against PubChem CID 31423 (pyrene), C16H10.
    assert smiles_to_iupac("c1cc2ccc3cccc4ccc(c1)c2c34") == "pyrene"


def test_acenaphthylene():
    # cross-checked against PubChem CID 9161 (acenaphthylene), C12H8; two
    # equally valid SMILES (aromatic and Kekulized bridge) both resolve.
    assert smiles_to_iupac("C1=Cc2cccc3cccc1c23") == "acenaphthylene"
    assert smiles_to_iupac("C1=CC2=CC=CC3=C2C(=C1)C=C3") == "acenaphthylene"


def test_substituted_pyrene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cc2ccc3cccc4ccc(c1)c2c34")


def test_substituted_acenaphthylene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2cccc3C=Cc1c23")


def test_unrelated_peri_fused_shape_raises():
    # a larger peri-fused hydrocarbon (verified via RDKit: some atom is
    # shared by three rings, but its formula/skeleton is neither pyrene
    # nor acenaphthylene) must still fall through to the ordinary
    # "not supported" path, not accidentally match one of them.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1cc2ccc3cc4ccc5ccc6cc1c1c2c3c4c5c61")
