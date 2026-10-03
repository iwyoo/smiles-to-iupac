import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified exactly.
        ("NCCC(=O)OC", "methyl 3-aminopropanoate"),
        ("NCC(=O)OCC", "ethyl 2-aminoethanoate"),
    ],
)
def test_smiles_to_iupac_ester_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_amine_on_alcohol_part():
    assert smiles_to_iupac("CC(=O)OCCN") == "2-aminoethyl ethanoate"


def test_two_amines():
    assert smiles_to_iupac("NCC(N)C(=O)OC") == "methyl 2,3-diaminopropanoate"


def test_two_esters_names_polyester():
    assert smiles_to_iupac("NCC(=O)OCOC(=O)C") == "methylene 2-aminoethanoate ethanoate"


def test_ring():
    assert smiles_to_iupac("NC1CCCCC1C(=O)OC") == "methyl 2-aminocyclohexane-1-carboxylate"


def test_unsaturated_chain():
    assert smiles_to_iupac("NCC=CC(=O)OC") == "methyl 4-aminobut-2-enoate"


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](C)C(=O)OC")


def test_hydroxyl_coexisting():
    assert smiles_to_iupac("NCC(O)C(=O)OC") == "methyl 3-amino-2-hydroxypropanoate"


def test_plain_ester_still_works():
    assert smiles_to_iupac("CC(=O)OC") == "methyl ethanoate"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"
