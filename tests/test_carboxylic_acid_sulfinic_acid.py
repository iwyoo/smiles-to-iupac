import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a)-style, no locant needed since
        # the -COOH carbon is the sole chain carbon). PubChem CID
        # 17889040 gives the retained name "sulfinoformic acid"; this
        # codebase always uses the systematic stem instead (see
        # `_carboxylic_acid_sulfonic_acid.py`'s identical precedent).
        ("OC(=O)S(=O)O", "sulfinomethanoic acid"),
        # PubChem CID 13101808 gives "2-sulfinoacetic acid"; same
        # systematic-stem substitution.
        ("OC(=O)CS(=O)O", "2-sulfinoethanoic acid"),
        # PubChem CID 3016737.
        ("OC(=O)CCS(=O)O", "3-sulfinopropanoic acid"),
        # PubChem CID 20563172.
        ("OC(=O)CCCS(=O)O", "4-sulfinobutanoic acid"),
    ],
)
def test_carboxylic_acid_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_carboxylic_acid_still_routes_normally():
    # No sulfinic acid present -- must still reach `_carboxylic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_plain_sulfinic_acid_still_routes_normally():
    # No carboxylic acid present -- must still reach `_sulfinic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)O") == "methanesulfinic acid"


def test_multiple_carboxylic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CC(S(=O)O)C(=O)O")


def test_other_heteroatom_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C(O)CS(=O)O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C=CCS(=O)O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC(S(=O)O)CC1")
