import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified exactly (CID 192802, CID 55290469).
        ("NCCC(N)=O", "3-aminopropanamide"),
        ("NCC(Cl)C(N)=O", "3-amino-2-chloropropanamide"),
    ],
)
def test_smiles_to_iupac_amide_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_alkyl_amide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(=O)CCN")


def test_two_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(N)C(N)=O")


def test_two_amides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)CC(N)=O")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCCC1C(N)=O")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC=CC(N)=O")


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](C)C(N)=O")


def test_hydroxyl_coexisting_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(O)C(N)=O")


def test_plain_amide_still_works():
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"
