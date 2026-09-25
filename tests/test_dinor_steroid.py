from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._dinor_steroid import _DINOR_LOOKUP


def _smiles_for(lo, hi, parent):
    label = f"{lo},{hi}"
    for smi, (lbl, name) in _DINOR_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no dinor-lookup entry for {label}-dinor-{parent}")


def test_lookup_spans_multiple_parent_families():
    names = {name for _, name in _DINOR_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_cross_ring_pair_androstane():
    assert smiles_to_iupac(_smiles_for(1, 12, "androstane")) == "1,12-dinor-androstane"


def test_both_methyls_excludes_as_gonane_collision():
    labels_names = {(lbl, name) for lbl, name in _DINOR_LOOKUP.values()}
    assert ("18,19", "androstane") not in labels_names


def test_ring_c_pair_cholestane():
    assert smiles_to_iupac(_smiles_for(11, 15, "cholestane")) == "11,15-dinor-cholestane"


def test_ring_pair_pregnane():
    assert smiles_to_iupac(_smiles_for(2, 6, "pregnane")) == "2,6-dinor-pregnane"


def test_gonane_pair():
    assert smiles_to_iupac(_smiles_for(1, 15, "gonane")) == "1,15-dinor-gonane"


def test_estrane_pair():
    assert smiles_to_iupac(_smiles_for(3, 11, "estrane")) == "3,11-dinor-estrane"


def test_18_19_dinor_androstane_resolves_to_gonane():
    from rdkit import Chem
    from rdkit.Chem import RWMol

    androstane = Chem.MolFromSmiles("CC12CCCC1C3CCC4CCCCC4(C3CC2)C")
    methyls = sorted([a.GetIdx() for a in androstane.GetAtoms() if a.GetDegree() == 1], reverse=True)
    rw = RWMol(androstane)
    for idx in methyls:
        rw.RemoveAtom(idx)
    m2 = rw.GetMol()
    Chem.SanitizeMol(m2)
    assert smiles_to_iupac(Chem.MolToSmiles(m2)) == "gonane"


def test_dinor_lookup_uses_lowest_locant_pair():
    for lbl, name in _DINOR_LOOKUP.values():
        a, b = (int(x) for x in lbl.split(","))
        assert a < b


def test_plain_ergostane_still_resolves():
    assert smiles_to_iupac("CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "ergostane"
