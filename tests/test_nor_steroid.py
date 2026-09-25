from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._nor_steroid import _NOR_LOOKUP


def _smiles_for(locant, parent):
    for smi, (loc, name) in _NOR_LOOKUP.items():
        if loc == locant and name == parent:
            return smi
    raise AssertionError(f"no nor-lookup entry for {locant}-nor-{parent}")


def test_lookup_spans_multiple_parent_families():
    names = {name for _, name in _NOR_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_4_nor_pregnane_blue_book_example_position():
    assert smiles_to_iupac(_smiles_for(4, "pregnane")) == "4-nor-pregnane"


def test_18_nor_androstane_ring_methyl():
    assert smiles_to_iupac(_smiles_for(18, "androstane")) == "18-nor-androstane"


def test_11_nor_cholestane_ring_position():
    assert smiles_to_iupac(_smiles_for(11, "cholestane")) == "11-nor-cholestane"


def test_6_nor_ergostane_ring_position():
    assert smiles_to_iupac(_smiles_for(6, "ergostane")) == "6-nor-ergostane"


def test_15_nor_cholane_ring_position():
    assert smiles_to_iupac(_smiles_for(15, "cholane")) == "15-nor-cholane"


def test_17_nor_gonane_ring_position():
    assert smiles_to_iupac(_smiles_for(17, "gonane")) == "17-nor-gonane"


def test_3_nor_estrane_ring_position():
    assert smiles_to_iupac(_smiles_for(3, "estrane")) == "3-nor-estrane"


def test_1_nor_androstane_ring_position():
    assert smiles_to_iupac(_smiles_for(1, "androstane")) == "1-nor-androstane"


def test_19_nor_androstane_collides_with_estrane_not_claimed_as_nor():
    androstane = Chem.MolFromSmiles("CC12CCCC1C3CCC4CCCCC4(C3CC2)C")
    rw = Chem.RWMol(androstane)
    rw.RemoveAtom(18)
    Chem.SanitizeMol(rw)
    assert smiles_to_iupac(Chem.MolToSmiles(rw)) == "estrane"


def test_18_nor_estrane_collides_with_gonane_not_claimed_as_nor():
    estrane = Chem.MolFromSmiles("CC12CCCC1C1CCC3CCCCC3C1CC2")
    methyl = [a.GetIdx() for a in estrane.GetAtoms() if a.GetDegree() == 1][0]
    rw = Chem.RWMol(estrane)
    rw.RemoveAtom(methyl)
    Chem.SanitizeMol(rw)
    assert smiles_to_iupac(Chem.MolToSmiles(rw)) == "gonane"


def test_ring_fusion_atom_not_eligible_for_nor():
    names = {(loc, name) for loc, name in _NOR_LOOKUP.values()}
    for parent in ("gonane", "androstane", "estrane", "pregnane", "cholane", "cholestane", "ergostane"):
        for fusion_locant in (5, 8, 9, 10, 13, 14):
            assert (fusion_locant, parent) not in names


def test_plain_androstane_still_resolves():
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"
