import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a), no locant on the -SO3H).
        # PubChem CID 54528747.
        ("NS(=O)(=O)CS(=O)(=O)O", "sulfamoylmethanesulfonic acid"),
        # Two-carbon chain: not the P-14.3.4.2(b) no-locant case, since a
        # sulfamoyl substituent is present in addition to the -SO3H
        # suffix (mirrors `_sulfonic_acid_sulfinic_acid.py`'s identical
        # -1- citation for the same chain shape; PubChem's own computed
        # name for this CID, 19002484, omits the locant, but PubChem's
        # algorithmic names are not always strict PINs).
        ("NS(=O)(=O)CCS(=O)(=O)O", "2-sulfamoylethanesulfonic acid"),
        # PubChem CID 119095391.
        ("NS(=O)(=O)CCCS(=O)(=O)O", "3-sulfamoylpropane-1-sulfonic acid"),
        # Halogen coexisting alongside the sulfonic acid/sulfonamide
        # pair. PubChem CID 19913272.
        ("NS(=O)(=O)CC(Cl)S(=O)(=O)O", "1-chloro-2-sulfamoylethanesulfonic acid"),
    ],
)
def test_sulfonic_acid_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_sulfonic_acid_still_routes_normally():
    # No sulfonamide present -- must still reach `_sulfonic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_plain_sulfonamide_still_routes_normally():
    # No sulfonic acid present -- must still reach `_sulfonamide.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"


def test_multiple_sulfonic_acids():
    assert smiles_to_iupac("OS(=O)(=O)CCS(=O)(=O)O") == "ethane-1,2-disulfonic acid"


def test_n_substituted_sulfonamide():
    assert smiles_to_iupac("CNS(=O)(=O)CCS(=O)(=O)O") == "2-(methylsulfamoyl)ethanesulfonic acid"


def test_other_heteroatom():
    assert smiles_to_iupac("OCC(S(=O)(=O)N)S(=O)(=O)O") == "2-hydroxy-1-sulfoethanesulfonamide"


def test_unsaturated_chain():
    assert smiles_to_iupac("C=CC(S(=O)(=O)N)S(=O)(=O)O") == "1-sulfoprop-2-ene-1-sulfonamide"


def test_ring():
    assert smiles_to_iupac("C1CCC(S(=O)(=O)N)(CC1)S(=O)(=O)O") == "1-sulfamoylcyclohexane-1-sulfonic acid"
