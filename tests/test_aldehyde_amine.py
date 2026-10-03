import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID for this exact structure gives '2-aminoacetaldehyde'
        # (the retained 'acetaldehyde' stem); this project's own
        # `_aldehyde.py` always uses the systematic 'ethanal' stem instead
        # (see that module's own tests), so this module follows that same
        # pre-existing convention.
        ("NCC=O", "2-aminoethanal"),
        # 3-aminopropanal: PubChem-verified exactly (no stem divergence for
        # a 3+ carbon chain).
        ("NCCC=O", "3-aminopropanal"),
        # A halogen substituent coexists with both the aldehyde and the
        # amine.
        ("NC(Cl)C=O", "2-amino-2-chloroethanal"),
    ],
)
def test_smiles_to_iupac_aldehyde_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_secondary_amine():
    assert smiles_to_iupac("CNCC=O") == "2-(methylamino)ethanal"


def test_two_amines():
    assert smiles_to_iupac("NC(N)C=O") == "2,2-diaminoethanal"


def test_two_aldehydes():
    assert smiles_to_iupac("O=CC(N)C=O") == "2-aminopropanedial"


def test_ring():
    assert smiles_to_iupac("NC1CCC(C=O)CC1") == "4-aminocyclohexane-1-carbaldehyde"


def test_unsaturated_chain():
    assert smiles_to_iupac("NCC=CC=O") == "4-aminobut-2-enal"


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](N)C=O")


def test_plain_aldehyde_still_works():
    assert smiles_to_iupac("CCC=O") == "propanal"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"


def test_hydroxyl_coexisting():
    assert smiles_to_iupac("NC(CO)C=O") == "2-amino-3-hydroxypropanal"
