import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a), no locant on the 'sulfamoyl'
        # prefix). PubChem CID 19097208 gives the retained name
        # "sulfamoylformic acid"; this codebase always uses the
        # systematic stem instead (see
        # `_carboxylic_acid_sulfinic_acid.py`'s identical precedent).
        ("OC(=O)S(=O)(=O)N", "sulfamoylmethanoic acid"),
        # PubChem CID 11083961 gives "2-sulfamoylacetic acid"; same
        # systematic-stem substitution.
        ("OC(=O)CS(=O)(=O)N", "2-sulfamoylethanoic acid"),
        # PubChem CID 21718767.
        ("OC(=O)CCS(=O)(=O)N", "3-sulfamoylpropanoic acid"),
        # PubChem CID 2773320.
        ("OC(=O)CCCS(=O)(=O)N", "4-sulfamoylbutanoic acid"),
    ],
)
def test_carboxylic_acid_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_carboxylic_acid_still_routes_normally():
    # No sulfonamide present -- must still reach `_carboxylic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_plain_sulfonamide_still_routes_normally():
    # No carboxylic acid present -- must still reach `_sulfonamide.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"


def test_multiple_carboxylic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CC(S(=O)(=O)N)C(=O)O")


def test_n_substituted_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CCS(=O)(=O)NC")


def test_other_heteroatom_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C(O)CS(=O)(=O)N")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C=CCS(=O)(=O)N")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC(S(=O)(=O)N)CC1")
