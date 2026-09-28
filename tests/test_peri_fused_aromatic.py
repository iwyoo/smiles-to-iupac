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


def test_fluoranthene():
    # cross-checked against PubChem CID 9154 (fluoranthene), C16H10 -- a
    # plain benzo ring ortho-fused onto acenaphthylene's own C1=C2 bond.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C3=CC=CC4=C3C2=CC=C4") == "fluoranthene"


def test_aceanthrylene():
    # cross-checked against PubChem CID 107781 (aceanthrylene), C16H10 --
    # a plain benzo ring ortho-fused onto acenaphthylene's ring, linear
    # arrangement (the C6C6C6C5 shape distinct from acephenanthrylene).
    assert smiles_to_iupac("C1=CC=C2C3=C4C(=CC=CC4=CC2=C1)C=C3") == "aceanthrylene"


def test_acephenanthrylene():
    # cross-checked against PubChem CID 9143 (acephenanthrylene), C16H10
    # -- the angular counterpart of aceanthrylene.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=C3C=CC4=C3C2=CC=C4") == "acephenanthrylene"


def test_coronene():
    # cross-checked against PubChem CID 9115 (coronene), C24H12; PIN per
    # tmp/bluebook/P2.txt line 2956.
    assert smiles_to_iupac("C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7") == "coronene"


def test_perylene():
    # cross-checked against PubChem CID 9142 (perylene), C20H12; PIN per
    # tmp/bluebook/P2.txt line 2962.
    assert smiles_to_iupac("C1=CC2=C3C(=C1)C4=CC=CC5=C4C(=CC=C5)C3=CC=C2") == "perylene"


def test_substituted_pyrene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cc2ccc3cccc4ccc(c1)c2c34")


def test_substituted_acenaphthylene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2cccc3C=Cc1c23")


def test_substituted_fluoranthene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)-c1cccc3cccc-2c13")


def test_unrelated_peri_fused_shape_raises():
    # a larger peri-fused hydrocarbon (verified via RDKit: some atom is
    # shared by three rings, but its formula/skeleton is neither pyrene
    # nor acenaphthylene) must still fall through to the ordinary
    # "not supported" path, not accidentally match one of them.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1cc2ccc3cc4ccc5ccc6cc1c1c2c3c4c5c61")
