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


def test_multiple_carboxylic_acids():
    assert smiles_to_iupac("OC(=O)CC(S(=O)(=O)N)C(=O)O") == "2-sulfamoylbutanedioic acid"


def test_n_substituted_sulfonamide():
    assert smiles_to_iupac("OC(=O)CCS(=O)(=O)NC") == "3-(methylsulfamoyl)propanoic acid"


def test_other_heteroatom():
    assert smiles_to_iupac("OC(=O)C(O)CS(=O)(=O)N") == "2-hydroxy-3-sulfamoylpropanoic acid"


def test_unsaturated_chain():
    assert smiles_to_iupac("OC(=O)C=CCS(=O)(=O)N") == "4-sulfamoylbut-2-enoic acid"


def test_ring():
    assert smiles_to_iupac("OC(=O)C1CCC(S(=O)(=O)N)CC1") == "4-sulfamoylcyclohexane-1-carboxylic acid"


def test_phenyl_chain_carboxylic_acid_sulfonamide():
    # A plain, unsubstituted benzene ring on the carboxylic acid/
    # sulfonamide chain, mirroring `_carboxylic_acid_sulfonic_acid.py`'s
    # own benzene-ring-substituent path (PR #362), cross-checked against
    # PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1CC(S(=O)(=O)N)C(=O)O") == "3-phenyl-2-sulfamoylpropanoic acid"  # CID 70062822


def test_phenyl_ring_with_second_substituent():
    assert smiles_to_iupac("OC(=O)c1ccccc1S(=O)(=O)N") == "2-sulfamoylbenzoic acid"


def test_phenyl_chain_carboxylic_acid_sulfonamide_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(S(=O)(=O)N)C(=O)O") == "3-(2-ethenylphenyl)-2-sulfamoylpropanoic acid"
