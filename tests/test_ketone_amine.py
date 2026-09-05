import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified exactly (no retained-name divergence for this
        # class, unlike the alcohol/thiol/aldehyde pairwise modules).
        ("NCC(C)=O", "1-aminopropan-2-one"),
        ("NCCC(C)=O", "4-aminobutan-2-one"),
        # A halogen substituent coexists with both the ketone and the
        # amine.
        ("NC(Cl)C(C)=O", "1-amino-1-chloropropan-2-one"),
    ],
)
def test_smiles_to_iupac_ketone_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_secondary_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNCC(C)=O")


def test_two_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(N)C(C)=O")


def test_two_ketones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(=O)CC(C)=O")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCC(=O)CC1")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC=CC(C)=O")


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](N)C(C)=O")


def test_plain_ketone_still_works():
    assert smiles_to_iupac("CC(C)=O") == "propan-2-one"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"


def test_hydroxyl_coexisting_raises():
    # A standalone hydroxyl alongside the ketone/amine pair is out of
    # scope for this first pairwise pilot on `_ketone.py`.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(=O)CO")
