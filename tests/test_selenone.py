import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Diethyl selenone (R=R'=ethyl), structure confirmed on PubChem
        # (CID 59266732, auto-generated non-PIN name '1-ethylselenonylethane').
        # PIN uses the acid-derived acyl prefix 'ethaneselenonyl' per
        # P-63.6, mirroring `_sulfone.py`'s '(ethanesulfonyl)ethane' exactly
        # with Se in place of S -- two-carbon symmetric case, no locant
        # (P-14.3.4.2(b)).
        ("CC[Se](=O)(=O)CC", "(ethaneselenonyl)ethane"),
        # Symmetric mononuclear parent (P-14.3.4.2(a), no locant): dimethyl
        # selenone, structure confirmed on PubChem (CID 176413, auto name
        # 'methylselenonylmethane').
        ("C[Se](=O)(=O)C", "(methaneselenonyl)methane"),
        # Asymmetric case: the longer chain (butane, R) is the parent, the
        # shorter (methyl, R') becomes the acid-derived acyl prefix.
        ("C[Se](=O)(=O)CCCC", "1-(methaneselenonyl)butane"),
    ],
)
def test_selenone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se](=O)(=O)C")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Se]1(=O)CCCCC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Se](=O)(=O)C")


def test_two_selenone_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Se](=O)(=O)C[Se](=O)(=O)C")
