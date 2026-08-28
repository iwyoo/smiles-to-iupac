import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Ethyl butyl selenoxide (R=butyl, R'=ethyl), structure confirmed on
        # PubChem (CID 13640793, whose own auto-generated non-PIN name is
        # '1-ethylseleninylbutane'). PIN uses the acid-derived acyl prefix
        # 'ethaneseleninyl' per P-63.6, mirroring `_sulfoxide.py`'s
        # '1-(ethanesulfinyl)butane' exactly with Se in place of S.
        ("CCCC[Se](=O)CC", "1-(ethaneseleninyl)butane"),
        # Symmetric case, mononuclear parent (P-14.3.4.2(a), no locant):
        # dimethyl selenoxide, structure confirmed on PubChem (CID 134574,
        # auto-generated non-PIN name 'methylseleninylmethane').
        ("C[Se](C)=O", "(methaneseleninyl)methane"),
    ],
)
def test_selenoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se](=O)C")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Se]1CCCCC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Se](=O)C")


def test_two_selenoxide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Se](=O)C[Se](=O)C")
