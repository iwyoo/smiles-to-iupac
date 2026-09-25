from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._homo_steroid import _HOMO_LOOKUP
from smiles_to_iupac._nor_steroid import _NOR_LOOKUP


def _smiles_for(label, parent):
    for smi, (lbl, name) in _HOMO_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no homo-lookup entry for {label}-homo-{parent}")


def test_lookup_spans_multiple_parent_families():
    names = {name for _, name in _HOMO_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_bond_connector_blue_book_worked_example():
    assert smiles_to_iupac(_smiles_for("13(17)a", "pregnane")) == "13(17)a-homo-pregnane"


def test_terminal_methyl_blue_book_worked_example():
    assert smiles_to_iupac(_smiles_for("19a", "pregnane")) == "19a-homo-pregnane"


def test_terminal_methyl_c18():
    assert smiles_to_iupac(_smiles_for("18a", "androstane")) == "18a-homo-androstane"


def test_ring_a_bond_androstane():
    assert smiles_to_iupac(_smiles_for("1(2)a", "androstane")) == "1(2)a-homo-androstane"


def test_ring_fusion_bond_cholestane():
    assert smiles_to_iupac(_smiles_for("8(9)a", "cholestane")) == "8(9)a-homo-cholestane"


def test_ring_c_bond_ergostane():
    assert smiles_to_iupac(_smiles_for("9(11)a", "ergostane")) == "9(11)a-homo-ergostane"


def test_ring_b_bond_cholane():
    assert smiles_to_iupac(_smiles_for("5(6)a", "cholane")) == "5(6)a-homo-cholane"


def test_ring_d_bond_estrane():
    assert smiles_to_iupac(_smiles_for("13(17)a", "estrane")) == "13(17)a-homo-estrane"


def test_gonane_lowest_locant_tie_break():
    assert smiles_to_iupac(_smiles_for("1(2)a", "gonane")) == "1(2)a-homo-gonane"


def test_plain_pregnane_still_resolves():
    assert smiles_to_iupac("CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "pregnane"


def test_nor_androstane_still_resolves():
    nor_smiles = next(s for s, (loc, name) in _NOR_LOOKUP.items() if loc == 3 and name == "androstane")
    assert smiles_to_iupac(nor_smiles) == "3-nor-androstane"
