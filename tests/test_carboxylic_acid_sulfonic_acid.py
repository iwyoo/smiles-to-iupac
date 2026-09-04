import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 31257 gives the retained-name "2-sulfoacetic acid";
        # this codebase always uses the systematic '...oic acid' stem
        # instead of a retained name once any substituent is present
        # (see `_carboxylic_acid.py`'s own established convention, e.g.
        # 'ClCC(=O)O' -> '2-chloroethanoic acid' not '2-chloroacetic
        # acid'), so 'ethanoic' replaces 'acetic' here too.
        ("OC(=O)CS(=O)(=O)O", "2-sulfoethanoic acid"),
        # PubChem CID 409694.
        ("OC(=O)CCS(=O)(=O)O", "3-sulfopropanoic acid"),
        # PubChem CID 21206169.
        ("OC(=O)CCCS(=O)(=O)O", "4-sulfobutanoic acid"),
        # Halogen coexisting alongside the pair. PubChem CID 4308193 gives
        # "2-chloro-2-sulfoacetic acid"; same systematic-stem substitution
        # as above.
        ("OC(=O)C(Cl)S(=O)(=O)O", "2-chloro-2-sulfoethanoic acid"),
        # Multiple sulfonic acids (multiplying prefix). PubChem CID
        # 89390889.
        ("OC(=O)CC(S(=O)(=O)O)CS(=O)(=O)O", "3,4-disulfobutanoic acid"),
    ],
)
def test_carboxylic_acid_sulfonic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_carboxylic_acid_still_routes_normally():
    # No sulfonic acid present -- must still reach `_carboxylic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_plain_sulfonic_acid_still_routes_normally():
    # No carboxylic acid present -- must still reach `_sulfonic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_multiple_carboxylic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CC(S(=O)(=O)O)C(=O)O")


def test_other_heteroatom_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C(O)CS(=O)(=O)O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C=CCS(=O)(=O)O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC(S(=O)(=O)O)CC1")
