from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._cyclo_steroid import _CYCLO_LOOKUP


def _smiles_for(lo, hi, parent):
    label = f"{lo},{hi}"
    for smi, (lbl, name) in _CYCLO_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no cyclo-lookup entry for {label}-cyclo-{parent}")


def test_lookup_spans_multiple_parent_families():
    names = {name for _, name in _CYCLO_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_blue_book_worked_example_shape():
    assert smiles_to_iupac(_smiles_for(3, 5, "pregnane")) == "3,5-cyclo-pregnane"


def test_fusion_atom_terminus_androstane():
    assert smiles_to_iupac(_smiles_for(3, 5, "androstane")) == "3,5-cyclo-androstane"


def test_ring_c_pair_cholestane():
    assert smiles_to_iupac(_smiles_for(9, 12, "cholestane")) == "9,12-cyclo-cholestane"


def test_cross_ring_pair_ergostane():
    assert smiles_to_iupac(_smiles_for(1, 6, "ergostane")) == "1,6-cyclo-ergostane"


def test_gonane_pair():
    assert smiles_to_iupac(_smiles_for(2, 7, "gonane")) == "2,7-cyclo-gonane"


def test_estrane_pair():
    assert smiles_to_iupac(_smiles_for(3, 6, "estrane")) == "3,6-cyclo-estrane"


def test_no_bond_pair_excluded():
    labels_names = set(_CYCLO_LOOKUP.values())
    assert ("1,2", "androstane") not in labels_names


def test_methyl_locants_excluded():
    for label, _ in _CYCLO_LOOKUP.values():
        a, b = (int(x) for x in label.split(","))
        assert a <= 17 and b <= 17


def test_plain_cholane_still_resolves():
    assert smiles_to_iupac("CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholane"


def test_seco_androstane_still_resolves():
    from smiles_to_iupac._seco_steroid import _SECO_LOOKUP

    seco_smiles = next(s for s, (lbl, name) in _SECO_LOOKUP.items() if lbl == "9,10" and name == "androstane")
    assert smiles_to_iupac(seco_smiles) == "9,10-seco-androstane"
