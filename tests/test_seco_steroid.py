from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._seco_steroid import _SECO_LOOKUP


def _smiles_for(lo, hi, parent):
    label = f"{lo},{hi}"
    for smi, (lbl, name) in _SECO_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no seco-lookup entry for {label}-seco-{parent}")


def test_lookup_spans_multiple_parent_families():
    names = {name for _, name in _SECO_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_plain_ring_bond_cleavage():
    assert smiles_to_iupac(_smiles_for(2, 3, "androstane")) == "2,3-seco-androstane"


def test_vitamin_d_anchor_ring_fusion_bond_androstane():
    assert smiles_to_iupac(_smiles_for(9, 10, "androstane")) == "9,10-seco-androstane"


def test_vitamin_d_anchor_ring_fusion_bond_cholestane():
    assert smiles_to_iupac(_smiles_for(9, 10, "cholestane")) == "9,10-seco-cholestane"


def test_ring_c_bond_ergostane():
    assert smiles_to_iupac(_smiles_for(11, 12, "ergostane")) == "11,12-seco-ergostane"


def test_ring_d_bond_pregnane():
    assert smiles_to_iupac(_smiles_for(13, 17, "pregnane")) == "13,17-seco-pregnane"


def test_ring_b_bond_cholane():
    assert smiles_to_iupac(_smiles_for(6, 7, "cholane")) == "6,7-seco-cholane"


def test_gonane_ring_a_bond():
    assert smiles_to_iupac(_smiles_for(3, 4, "gonane")) == "3,4-seco-gonane"


def test_estrane_ring_fusion_bond():
    assert smiles_to_iupac(_smiles_for(8, 14, "estrane")) == "8,14-seco-estrane"


def test_plain_cholestane_still_resolves():
    assert smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholestane"


def test_homo_androstane_still_resolves():
    from smiles_to_iupac._homo_steroid import _HOMO_LOOKUP

    homo_smiles = next(s for s, (lbl, name) in _HOMO_LOOKUP.items() if lbl == "19a" and name == "androstane")
    assert smiles_to_iupac(homo_smiles) == "19a-homo-androstane"
