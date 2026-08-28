import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-63.6's own worked example for diethyl sulfone
        # (R=R'=ethyl): '(ethanesulfonyl)ethane (PIN)' -- note this is the
        # two-carbon symmetric case with NO locant despite the enclosed
        # acyl-derived prefix (P-14.3.4.2(b), see `_sulfoxide.py`'s module
        # docstring for why).
        ("CCS(=O)(=O)CC", "(ethanesulfonyl)ethane"),
        # Symmetric mononuclear parent (P-14.3.4.2(a), no locant): dimethyl
        # sulfone. Differs from PubChem's non-PIN 'methylsulfonylmethane'
        # for the same reason `_sulfoxide.py`'s DMSO test does.
        ("CS(=O)(=O)C", "(methanesulfonyl)methane"),
        # Asymmetric case: the longer chain (butane, R) is the parent, the
        # shorter (methyl, R') becomes the acid-derived acyl prefix.
        ("CS(=O)(=O)CCCC", "1-(methanesulfonyl)butane"),
    ],
)
def test_sulfone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)S(=O)(=O)C")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S1(=O)CCCCC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CS(=O)(=O)C")


def test_two_sulfone_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)(=O)CS(=O)(=O)C")
