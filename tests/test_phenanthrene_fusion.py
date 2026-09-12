import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_benzo_c_phenanthrene():
    # PubChem CID 9136.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C=C3") == "benzo[c]phenanthrene"


def test_chrysene():
    # Letter 'a' -- chrysene (CID 9171) -- is a retained name, cited as
    # such rather than the systematic 'benzo[a]phenanthrene'.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43") == "chrysene"


def test_letter_b_defers_to_anthracene_module():
    # Letter 'b' on phenanthrene is the exact same compound as letter 'a'
    # on anthracene (benzo[a]anthracene, CID 5954) -- must be excluded
    # here so `_anthracene_fusion.py` (the senior base component) claims
    # it instead, not duplicated or shadowed.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32") == "benzo[a]anthracene"


def test_letter_l_defers_to_triphenylene():
    # Letter 'l' (the K-region, 9,10-bond) is triphenylene (CID 9170),
    # already a retained name via `_branched_fused_aromatic.py`.
    assert smiles_to_iupac("c1ccc2c(c1)c1ccccc1c1ccccc21") == "triphenylene"


def test_substituted_chrysene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)ccc1c2ccc2ccccc21")
