import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Symmetric case, mononuclear parent (P-14.3.4.2(a), no locant):
        # dimethyl telluroxide, structure confirmed on PubChem (CID
        # 14009075, auto-generated non-PIN name 'methyltellurinylmethane').
        ("C[Te](C)=O", "(methanetellurinyl)methane"),
        # Ethyl butyl telluroxide (R=butyl, R'=ethyl) -- no PubChem-registered
        # structure exists for this asymmetric case (CID 0), so this is a
        # reviewed, not directly verified, mechanical extension of the
        # identical O/S/Se pattern (see module docstring).
        ("CCCC[Te](=O)CC", "1-(ethanetellurinyl)butane"),
    ],
)
def test_telluroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te](=O)C")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Te]1CCCCC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Te](=O)C")


def test_two_telluroxide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Te](=O)C[Te](=O)C")
