import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_benzo_b_chrysene():
    # PubChem CID 9163 (a real registered compound, no trivial/retained
    # name of its own -- only a von Baeyer-style systematic name).
    assert smiles_to_iupac("c1ccc2cc3c(ccc4c5ccccc5ccc34)cc2c1") == "benzo[b]chrysene"


def test_benzo_c_chrysene():
    # PubChem CID 9135.
    assert smiles_to_iupac("c1ccc2c(c1)ccc1c2ccc2ccc3ccccc3c21") == "benzo[c]chrysene"


def test_benzo_g_chrysene():
    # PubChem CID 9140.
    assert smiles_to_iupac("c1ccc2c(c1)ccc1c3ccccc3c3ccccc3c21") == "benzo[g]chrysene"


def test_letter_a_defers_to_picene():
    # Letter 'a' on chrysene is picene (PubChem CID 9162), a retained
    # name -- excluded here, but the shape isn't supported anywhere else
    # yet either (picene itself has no dedicated recognizer), so it
    # currently still raises.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC=CC=C54")


def test_substituted_chrysene_fusion_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2cc3c(ccc4c5ccccc5ccc34)cc2c1")
