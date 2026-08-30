import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-64.6.1 worked example: "propane-2-thione (PIN) (not
        # thioacetone)". PubChem structure match confirms.
        ("CC(=S)C", "propane-2-thione"),
        # Blue Book P-64.6.1 worked example: "butane-2-thione (PIN)".
        ("CCC(=S)C", "butane-2-thione"),
        ("CC(=S)CC", "butane-2-thione"),
        # PubChem structure match: "hexane-3-thione".
        ("CCCC(=S)CC", "hexane-3-thione"),
        # Blue Book P-64.6.1 worked example: "pentane-2,4-dithione (PIN)".
        ("CC(=S)CC(=S)C", "pentane-2,4-dithione"),
        # PubChem structure match: "cyclohexanethione".
        ("C1CCC(=S)CC1", "cyclohexanethione"),
        # Halogen coexistence, PubChem structure match:
        # "1-chloropropane-2-thione".
        ("ClCC(=S)C", "1-chloropropane-2-thione"),
        # Unsaturated chain, PubChem structure match: "pent-4-ene-2-thione".
        ("C=CCC(=S)C", "pent-4-ene-2-thione"),
    ],
)
def test_thione_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_thial_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=S")


def test_aromatic_thione_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc(cc1)C(=S)C")


def test_polycyclic_thione_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CCC2(CCCCC2)CC1")


def test_thione_with_hydroxyl_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=S)C")


def test_thione_with_ketone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC(=S)C")


def test_thiol_not_confused_with_thione():
    assert smiles_to_iupac("CS") == "methanethiol"


def test_sulfide_not_confused_with_thione():
    assert smiles_to_iupac("CSC") == "methylsulfanylmethane"


def test_disulfide_not_confused_with_thione():
    assert smiles_to_iupac("CSSC") == "(methyldisulfanyl)methane"


def test_ketone_not_confused_with_thione():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"
