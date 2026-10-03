import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Ethanolamine: PubChem CID 700's own IUPACName is
        # '2-aminoethanol' (the retained 'ethanol' stem kept even when
        # substituted); this project's own `_alcohol.py` already diverges
        # from that convention for a substituted ethanol (e.g. its own
        # '2-ethoxyethan-1-ol' test), so this module follows that same
        # pre-existing convention instead.
        ("NCCO", "2-aminoethanol"),
        # 3-aminopropan-1-ol: PubChem-verified exactly (CID 7269).
        ("NCCCO", "3-aminopropan-1-ol"),
        # 2-aminopropan-1-ol: PubChem-verified exactly (the -OH gets the
        # lower locant, P-44).
        ("CC(N)CO", "2-aminopropan-1-ol"),
        # A second hydroxyl coexists with the amine -- PubChem-verified
        # exactly ('3-aminopropane-1,2-diol').
        ("NCC(O)CO", "3-aminopropane-1,2-diol"),
        # A halogen substituent coexists with both the hydroxyl and the
        # amine.
        ("NC(Cl)CO", "2-amino-2-chloroethanol"),
    ],
)
def test_smiles_to_iupac_alcohol_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_secondary_amine():
    assert smiles_to_iupac("CNCCO") == "2-(methylamino)ethanol"


def test_two_amines():
    assert smiles_to_iupac("NC(N)CO") == "2,2-diaminoethanol"


def test_ring():
    assert smiles_to_iupac("NC1CCC(O)CC1") == "4-aminocyclohexan-1-ol"


def test_unsaturated_chain():
    assert smiles_to_iupac("NCC=CCO") == "4-aminobut-2-en-1-ol"


def test_specified_stereocenter():
    assert smiles_to_iupac("C[C@H](N)CO") == "(2S)-2-aminopropan-1-ol"


def test_plain_alcohol_still_works():
    assert smiles_to_iupac("CCCO") == "propan-1-ol"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"


def test_ether_coexisting():
    assert smiles_to_iupac("NCC(O)COCC") == "1-amino-3-ethoxypropan-2-ol"
