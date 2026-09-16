import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-41 Table 4.1: hydroperoxide (class 18) outranks amine (class
        # 19) -- Blue Book primary source, followed directly here even
        # though PubChem's own auto-namer returns the opposite seniority
        # conclusion for this pair ("2-hydroperoxyethanamine" for
        # 'NCCOO'); see _hydroperoxide_amine.py's module docstring.
        ("NCCOO", "2-aminoethaneperoxol"),
        ("NCCCOO", "3-aminopropane-1-peroxol"),
    ],
)
def test_smiles_to_iupac_hydroperoxide_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(N)COO")


def test_two_hydroperoxides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(OO)COO")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCCC1COO")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC=CCOO")


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](C)COO")


def test_hydroxyl_coexisting_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(O)COO")


def test_plain_hydroperoxide_still_works():
    assert smiles_to_iupac("CCOO") == "ethaneperoxol"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"
