"""Diorganylphosphanyl groups (R2P-, P-68.3) as named monovalent prefixes.
`contract_phosphanyl_groups` swaps each such group for an iodine placeholder
carrying PREFIX_PROP, which `halogen_substituents` reads, so every halogen-
aware module cites e.g. '2-(diphenylphosphanyl)phenyl' without learning P.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._substituents import format_mononuclear_prefixes, name_branch

PREFIX_PROP = "_named_prefix"


def _branch_atoms(graph, root, parent):
    seen = {root}
    stack = [root]
    while stack:
        for v in graph[stack.pop()]:
            if v != parent and v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def phosphanyl_name(mol, graph, phosphorus, parent):
    """'diphenylphosphanyl' for the R2P- group on `phosphorus` reached from `parent`."""
    atom = mol.GetAtomWithIdx(phosphorus)
    if atom.GetFormalCharge() or atom.GetTotalNumHs() or atom.IsInRing():
        raise UnsupportedStructure("only a neutral acyclic R2P- group is supported as a prefix")
    roots = [n for n in graph[phosphorus] if n != parent]
    if len(roots) != 2 or any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 for n in roots):
        raise UnsupportedStructure("only a diorganylphosphanyl group is supported as a prefix")
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [name_branch(graph, n, phosphorus, {}, aromatic, mol=mol) for n in roots]
    return format_mononuclear_prefixes(entries) + "phosphanyl"


def contract_phosphanyl_groups(mol):
    """`mol` with each R2P- group collapsed to a placeholder, or None."""
    graph = adjacency(mol)
    groups = []
    for atom in mol.GetAtoms():
        p = atom.GetIdx()
        if atom.GetAtomicNum() != 15 or atom.GetDegree() != 3:
            continue
        sides = {n: _branch_atoms(graph, n, p) for n in graph[p]}
        if any(p in side for side in sides.values()):
            continue
        ordered = sorted(sides, key=lambda n: len(sides[n]), reverse=True)
        if len(sides[ordered[0]]) == len(sides[ordered[1]]):
            continue
        try:
            name = phosphanyl_name(mol, graph, p, ordered[0])
        except UnsupportedStructure:
            continue
        groups.append((p, ordered[0], {p} | sides[ordered[1]] | sides[ordered[2]], name))
    if not groups:
        return None
    rw = Chem.RWMol(mol)
    drop: set[int] = set()
    for p, parent, atoms, name in groups:
        placeholder = rw.AddAtom(Chem.Atom(53))
        rw.AddBond(parent, placeholder, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(placeholder).SetProp(PREFIX_PROP, f"({name})")
        drop |= atoms
    for idx in sorted(drop, reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out
