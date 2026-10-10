"""The nor/homo/seco search (P-101.3.7) recovers skeletons built by known operations, within the modifications allowed."""

import pytest
from rdkit import Chem

from smiles_to_iupac._np import _cyclomatic, _plausible, _view
from smiles_to_iupac._np_core import get_parent
from smiles_to_iupac._np_diff import WORK_LIMIT, WORK_PER_PARENT, read_operations
from smiles_to_iupac._np_match import View, embeddings
from smiles_to_iupac._np_region import regions
from smiles_to_iupac._np_skel import variants


def _molecule(skel):
    labels = sorted(skel.adj, key=str)
    index = {label: i for i, label in enumerate(labels)}
    mol = Chem.RWMol()
    table = Chem.GetPeriodicTable()
    for label in labels:
        mol.AddAtom(Chem.Atom(table.GetAtomicNumber(skel.elem[label])))
    for a, neighbours in skel.adj.items():
        for b in neighbours:
            if index[a] < index[b]:
                mol.AddBond(index[a], index[b], Chem.BondType.SINGLE)
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))


def _ops(skel):
    return tuple(op[:-1] if op[0] == "homo" else op for op in skel.ops)


@pytest.mark.parametrize("name, ops", [
    ("androstane", (("nor", "4"),)),
    ("morphinan", (("nor", "8"),)),
    ("yohimban", (("nor", "6"),)),
    ("2,2′-neolignane", (("nor", "6′"),)),
    ("podocarpane", (("nor", "15"), ("nor", "16"))),
    ("abietane", (("nor", "16"), ("nor", "17"))),
    ("3,7′-neolignane", (("homo", "bond", ("3", "7′")),)),
    ("cinchonan", (("homo", "atomic", ("6", "1")), ("seco", "4", "7"))),
    ("androstane", (("nor", "17"), ("seco", "13", "14"))),
    ("morphinan", (("nor", "4"), ("seco", "9", "10"))),
    ("yohimban", (("nor", "19"), ("seco", "3", "14"))),
    ("evonine", (("seco", "11", "12"),)),
], ids=lambda value: str(value) if isinstance(value, str) else "+".join(op[0] for op in value))
def test_search_recovers_the_operations_that_built_the_skeleton(name, ops):
    parent = get_parent(name)
    skel = next(s for s in variants(parent, len(ops)) if _ops(s) == ops)
    view = View(_molecule(skel))
    view.cyclomatic = _cyclomatic(view.adj)
    view.work = 10**9
    found = read_operations(parent, view, len(ops))
    assert found and any(embeddings(s, view, limit=1) for s in found)


def _skeleton(name, ops, replaced=()):
    skel = next(s for s in variants(get_parent(name), len(ops)) if _ops(s) == ops)
    for label in replaced:
        skel.elem[label] = "N"
    return skel


@pytest.mark.parametrize("replaced, expected", [(("1", "7", "12"), True), (("1", "7", "12", "16"), False)])
def test_search_leaves_no_more_replacements_than_modifications_remain(replaced, expected):
    # P-101.3.7.1: a nor and its replacements are all modifications, at most four
    parent = get_parent("androstane")
    view = View(_molecule(_skeleton("androstane", (("nor", "4"),), replaced)))
    view.cyclomatic = _cyclomatic(view.adj)
    view.work = 10**9
    found = read_operations(parent, view, 1)
    assert bool(found) is expected


def test_no_seco_is_offered_to_a_parent_with_two_rings():
    parent = get_parent("2,2′-neolignane")
    skel = next(s for s in variants(parent, 1) if s.ops[0][0] == "seco")
    view = View(_molecule(skel))
    view.cyclomatic = _cyclomatic(view.adj)
    view.work = 10**9
    assert not any(op[0] == "seco" for found in read_operations(parent, view, 1) for op in found.ops)


def test_search_of_a_parent_the_molecule_lacks_does_not_spend_the_work_budget():
    mol = Chem.MolFromSmiles("O=C(NCc1ccc2c(c1)OCO2)[C@@H]1CC(=O)N(c2ccc3c(c2)OCCO3)C1")
    view = _view(mol)
    assert read_operations(get_parent("18-oxayohimban"), view, 3) == []
    assert WORK_LIMIT - view.work < WORK_PER_PARENT


def test_regions_bridge_only_through_atoms_that_can_be_replaced():
    view = View(Chem.MolFromSmiles("CC(C)OC(C)C"))
    assert sorted(r.size for r in regions(view, {"C"}, 0)) == [3, 3]
    assert max(r.size for r in regions(view, {"C"}, 1)) == 7
    assert max(r.size for r in regions(View(Chem.MolFromSmiles("CC(C)C(C)O")), {"C"}, 1)) == 5


@pytest.mark.parametrize("replaced, expected", [(4, True), (5, False)])
def test_parent_needs_a_region_with_its_atoms_after_replacements(replaced, expected):
    skel = _skeleton("androstane", (), ("1", "2", "6", "7", "11")[:replaced])
    assert _plausible("androstane", _view(_molecule(skel)), 0) is expected


@pytest.mark.parametrize("stereo, expected", [(False, False), (True, True)])
def test_two_operations_need_a_stereo_element_in_the_region(stereo, expected):
    # a skeleton changed by two operations is named only when its atoms hold a specified configuration
    mol = _molecule(_skeleton("androstane", (("nor", "4"), ("nor", "17"))))
    if stereo:
        junction = next(a for a in mol.GetAtoms() if a.GetDegree() == 3)
        junction.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    assert _plausible("androstane", _view(mol), 2) is expected
