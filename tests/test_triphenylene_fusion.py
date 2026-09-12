import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_benzo_b_triphenylene():
    # PubChem CID 9164 -- PubChem's own autoname agrees independently.
    assert smiles_to_iupac("c1ccc2cc3c4ccccc4c4ccccc4c3cc2c1") == "benzo[b]triphenylene"


def test_letter_a_defers_to_chrysene_module():
    # Letter 'a' on triphenylene is the exact same compound as letter 'g'
    # on chrysene (benzo[g]chrysene, CID 9140) -- excluded here so
    # `_chrysene_fusion.py` claims it instead. Which base is actually
    # P-25.3.2.4-senior for this shape isn't fully verified (see
    # `_triphenylene_fusion.py`'s own module docstring); this pins down
    # current behavior, not a confirmed-correct name.
    assert smiles_to_iupac("c1ccc2c(c1)c1ccc3ccccc3c1c1ccccc21") == "benzo[g]chrysene"


def test_substituted_triphenylene_fusion_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2cc3c4ccccc4c4ccccc4c3cc2c1")
