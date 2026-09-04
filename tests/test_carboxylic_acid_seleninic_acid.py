import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 428823 gives the retained-name "2-seleninoacetic
        # acid"; this codebase always uses the systematic '...oic acid'
        # stem instead (see `_carboxylic_acid_sulfinic_acid.py`'s
        # identical precedent), so 'ethanoic' replaces 'acetic'.
        ("OC(=O)C[Se](=O)O", "2-seleninoethanoic acid"),
        # PubChem CID 99569.
        ("OC(=O)CC[Se](=O)O", "3-seleninopropanoic acid"),
        # PubChem CID 85843948.
        ("OC(=O)CCC[Se](=O)O", "4-seleninobutanoic acid"),
    ],
)
def test_carboxylic_acid_seleninic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mononuclear_parent_omits_locant():
    # No PubChem registration for this exact structure to check against,
    # but P-14.3.4.2(a) requires the sole substituent locant on a
    # one-carbon chain to be omitted -- same rule verified against
    # PubChem CID 18466245 ("sulfoformic acid") for the sulfonic-acid
    # sibling module (PR #315).
    assert smiles_to_iupac("OC(=O)[Se](=O)O") == "seleninomethanoic acid"


def test_plain_carboxylic_acid_still_routes_normally():
    # No seleninic acid present -- must still reach
    # `_carboxylic_acid.py`'s own path.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_plain_seleninic_acid_still_routes_normally():
    # No carboxylic acid present -- must still reach
    # `_seleninic_acid.py`'s own path.
    assert smiles_to_iupac("C[Se](=O)O") == "methaneseleninic acid"


def test_multiple_carboxylic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CC([Se](=O)O)C(=O)O")


def test_other_heteroatom_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C(O)C[Se](=O)O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C=CC[Se](=O)O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC([Se](=O)O)CC1")
