import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Symmetric mononuclear parent (P-14.3.4.2(a), no locant): dimethyl
        # tellurone, structure confirmed on PubChem (CID 59167377, auto name
        # 'methyltelluronylmethane').
        ("C[Te](=O)(=O)C", "(methanetelluronyl)methane"),
        # Diethyl tellurone -- no PubChem-registered structure exists for
        # this case (CID 0), so this is a reviewed, not directly verified,
        # mechanical extension of the identical O/S/Se pattern (see module
        # docstring). Two-carbon symmetric case, no locant (P-14.3.4.2(b)).
        ("CC[Te](=O)(=O)CC", "(ethanetelluronyl)ethane"),
        # Asymmetric case: the longer chain (butane, R) is the parent, the
        # shorter (methyl, R') becomes the acid-derived acyl prefix.
        ("C[Te](=O)(=O)CCCC", "1-(methanetelluronyl)butane"),
    ],
)
def test_tellurone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te](=O)(=O)C")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Te]1(=O)CCCCC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Te](=O)(=O)C")


def test_two_tellurone_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Te](=O)(=O)C[Te](=O)(=O)C")
