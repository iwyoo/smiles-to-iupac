"""A monocyclic metallacycle cited as a substituent group ('1-iodo-1-methyl-
...-1-platinacyclobutan-3-yl', P-69.4) of a more senior parent named by the
other modules; the ring is contracted to a placeholder carrying its name.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._metallacycle import FREE_VALENCE_PROP, _find_ring_metal, name_metallacycle
from ._phosphanyl_group import PREFIX_PROP
from ._prefix_groups import enclose


def _reach(graph, start, cut):
    seen, stack = {start}, [start]
    while stack:
        u = stack.pop()
        for v in graph[u]:
            if {u, v} == cut or v in seen:
                continue
            seen.add(v)
            stack.append(v)
    return seen


def name_metallacycle_as_group(mol):
    found = _find_ring_metal(mol)
    if found is None:
        raise UnsupportedStructure("no metallacycle here")
    metals, ring_atoms = found
    graph = adjacency(mol)
    ring = set(ring_atoms)
    from .core import _name_mol

    for a in ring_atoms:
        if a in metals:
            continue
        for x in graph[a]:
            if x in ring or mol.GetAtomWithIdx(x).GetAtomicNum() != 6:
                continue
            cut = {a, x}
            ring_side = _reach(graph, a, cut)
            if x in ring_side:
                continue
            group_mol = Chem.RWMol(mol)
            group_mol.GetAtomWithIdx(a).SetProp(FREE_VALENCE_PROP, "1")
            for idx in sorted(set(range(mol.GetNumAtoms())) - ring_side, reverse=True):
                group_mol.RemoveAtom(idx)
            group_mol = group_mol.GetMol()
            Chem.SanitizeMol(group_mol)
            try:
                group_name = name_metallacycle(group_mol)
            except UnsupportedStructure:
                continue
            rw = Chem.RWMol(mol)
            marker = rw.AddAtom(Chem.Atom(53))
            rw.AddBond(x, marker, Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(marker).SetProp(PREFIX_PROP, enclose(group_name))
            for idx in sorted(ring_side, reverse=True):
                rw.RemoveAtom(idx)
            rest = rw.GetMol()
            Chem.SanitizeMol(rest)
            try:
                name = _name_mol(rest)
            except UnsupportedStructure:
                continue
            if group_name in name and "iod" not in name.replace(group_name, ""):
                return name
    raise UnsupportedStructure("this metallacycle could not be cited as a substituent group")
