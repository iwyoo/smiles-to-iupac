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


def test_n_alkyl_amide():
    assert smiles_to_iupac("CNC(=O)CCN") == "3-amino-N-methylpropanamide"


def test_two_amines():
    assert smiles_to_iupac("NCC(N)C(N)=O") == "2,3-diaminopropanamide"


def test_two_amides():
    assert smiles_to_iupac("NC(=O)CC(N)=O") == "propanediamide"


def test_ring():
    assert smiles_to_iupac("NC1CCCCC1C(N)=O") == "2-aminocyclohexane-1-carboxamide"


def test_unsaturated_chain():
    assert smiles_to_iupac("NCC=CC(N)=O") == "4-aminobut-2-enamide"


def test_specified_stereocenter():
    assert smiles_to_iupac("N[C@@H](C)C(N)=O") == "(2S)-2-aminopropanamide"


def test_hydroxyl_coexisting():
    assert smiles_to_iupac("NCC(O)C(N)=O") == "3-amino-2-hydroxypropanamide"


def test_plain_amide_still_works():
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"
