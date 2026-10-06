"""The nor/homo/seco search (P-101.3.7) recovers skeletons that were built from a parent by known operations."""

import pytest
from rdkit import Chem

from smiles_to_iupac._np import _cyclomatic
from smiles_to_iupac._np_core import get_parent
from smiles_to_iupac._np_diff import read_operations
from smiles_to_iupac._np_match import View, embeddings
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
], ids=lambda value: str(value) if isinstance(value, str) else "+".join(op[0] for op in value))
def test_search_recovers_the_operations_that_built_the_skeleton(name, ops):
    parent = get_parent(name)
    skel = next(s for s in variants(parent, len(ops)) if _ops(s) == ops)
    view = View(_molecule(skel))
    view.cyclomatic = _cyclomatic(view.adj)
    view.work = 10**9
    found = read_operations(parent, view, len(ops))
    assert found and any(embeddings(s, view, limit=1) for s in found)
