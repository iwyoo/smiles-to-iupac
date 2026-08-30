import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem structure match: "propane-2-selone" (same structure as
        # the Blue Book's own "propane-2-thione (PIN)" for the sulfur case).
        ("CC(=[Se])C", "propane-2-selone"),
        # PubChem structure match: "pentane-2,4-diselone".
        ("CC(=[Se])CC(=[Se])C", "pentane-2,4-diselone"),
        # PubChem structure match: "cyclohexaneselone".
        ("C1CCC(=[Se])CC1", "cyclohexaneselone"),
        # PubChem structure match: "1-chloropropane-2-selone".
        ("ClCC(=[Se])C", "1-chloropropane-2-selone"),
    ],
)
def test_selone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_selenoaldehyde_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=[Se]")


def test_polycyclic_selone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCC2(CCCCC2)CC1")


def test_selone_with_hydroxyl_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=[Se])C")


def test_selenol_not_confused_with_selone():
    assert smiles_to_iupac("C[SeH]") == "methaneselenol"


def test_selenide_not_confused_with_selone():
    assert smiles_to_iupac("C[Se]C") == "methylselanylmethane"


def test_diselenide_not_confused_with_selone():
    assert smiles_to_iupac("C[Se][Se]C") == "(methyldiselanyl)methane"


def test_thione_not_confused_with_selone():
    assert smiles_to_iupac("CC(=S)C") == "propane-2-thione"


def test_ketone_not_confused_with_selone():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"
