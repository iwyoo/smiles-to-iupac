import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCCCC1CCC1CCCCC1", "(2-cyclohexylethyl)cyclohexane"),  # PubChem CID 76838
        ("C1CCCCC1CCCc1ccccc1", "(3-cyclohexylpropyl)benzene"),  # CID 561990
        ("c1ccccc1CCCc1ccccc1", "(3-phenylpropyl)benzene"),  # CID 14125
        ("C1CCCCC1CC1CCCCC1", "(cyclohexylmethyl)cyclohexane"),  # CID 76644
    ],
)
def test_disjoint_ring_pair_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_assembly_ylidene_still_raises():
    # Two disjoint rings joined by a C=C, not a plain single-bonded
    # bridge -- this module must not silently drop that double bond.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC1=C1CCCC1")


def test_three_disjoint_rings_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1CCC1CCCCC1CCC1CCCCC1")


def test_ring_with_extra_substituent_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CCCCC1CC1CCCCC1")
