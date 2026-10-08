"""Acyclic polyboranes with boron-boron bonds, and Lewis adducts of boranes (P-68.1.1.2.1, P-68.1.3, P-68.1.6).

The longest boron chain of n atoms is the multiplied parent hydride 'diborane', 'triborane', ... (other borons are
cited as boranyl branches); its hydrogen count, n + 2, follows in parentheses and survives an 'ene' or 'yne' ending.
A dative or charge-separated N(+)-B(-) pair is named as the adduct of its two neutral components (P-68.1.6.2)."""

from rdkit import Chem

from ._cited_group import cited_group
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, numerical_term, unsaturation_suffix
from ._hetero_prefixes import _has_senior_principal_group
from ._substituents import format_substituent_prefixes

_ACCEPTORS = {5, 13, 31}
_DONORS = {7, 8, 15, 16}


def _boron_paths(mol, borons):
    """Every longest path of the boron skeleton, in both directions, when the borons form one acyclic tree."""
    graph = {b: [n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetAtomicNum() == 5] for b in borons}
    if sum(len(v) for v in graph.values()) != 2 * (len(borons) - 1):
        return []
    reached, stack = set(), [borons[0]]
    while stack:
        node = stack.pop()
        if node not in reached:
            reached.add(node)
            stack.extend(graph[node])
    if len(reached) != len(borons):
        return []
    longest, paths = 0, []
    for start in borons:
        pending = [[start]]
        while pending:
            walk = pending.pop()
            onward = [n for n in graph[walk[-1]] if n not in walk]
            if onward:
                pending.extend(walk + [n] for n in onward)
            elif len(walk) > longest:
                longest, paths = len(walk), [walk]
            elif len(walk) == longest:
                paths.append(walk)
    return paths


def _valence(atom):
    return sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) + atom.GetTotalNumHs()


def _parent_name(count, ene, yne, substituted):
    """'diborane', 'triborene', 'tetrabor-1-ene': the ending of the chain hydride is modified (P-68.1.3, P-31.1.4)."""
    if not ene and not yne:
        return f"{numerical_term(count)}borane"
    root = f"{numerical_term(count)}bor"
    if count == 2 or (count == 3 and len(ene) + len(yne) == 1 and not substituted):
        return root + ("yne" if yne and not ene else "ene" if ene and not yne else unsaturation_suffix(ene, yne)[0])
    body, needs_a = unsaturation_suffix(ene, yne)
    return f"{root}{'a' if needs_a else ''}-{body}"


def _path_candidate(mol, graph, path):
    path_set = set(path)
    entries = []
    for position, boron in enumerate(path):
        for n in graph[boron]:
            if n in path_set:
                continue
            if mol.GetBondBetweenAtoms(boron, n).GetBondTypeAsDouble() != 1.0:
                return None
            entries.append((position + 1, *cited_group(mol, graph, n, boron)))
    ene, yne = [], []
    for i, (a, b) in enumerate(zip(path, path[1:])):
        order = int(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble())
        if order > 1:
            (ene if order == 2 else yne).append(i + 1)
    ordered = [locant for locant, name, _ in sorted(entries, key=lambda item: alpha_sort_key(item[1]))]
    return (ene, yne, sorted(locant for locant, _, _ in entries), ordered), entries, ene, yne


def polyborane_name(mol):
    borons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 5]
    if len(borons) < 2 or len(Chem.GetMolFrags(mol)) != 1 or _has_senior_principal_group(mol):
        return None
    if any(mol.GetAtomWithIdx(b).GetFormalCharge() or mol.GetAtomWithIdx(b).GetIsotope() for b in borons):
        return None
    if any(_valence(mol.GetAtomWithIdx(b)) != 3 for b in borons):
        return None
    if any(bond.GetBondTypeAsDouble() not in (1.0, 2.0, 3.0) for bond in mol.GetBonds()):
        return None
    paths = _boron_paths(mol, borons)
    if not paths:
        return None
    graph = adjacency(mol)
    best = None
    try:
        for path in paths:
            candidate = _path_candidate(mol, graph, path)
            if candidate is None:
                return None
            if best is None or candidate[0] < best[0]:
                best = (*candidate, len(path))
    except UnsupportedStructure:
        return None
    _, entries, ene, yne, count = best
    grouped = {}
    for locant, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    for info in grouped.values():
        info["locants"].sort()
    parent = _parent_name(count, ene, yne, bool(grouped))
    return f"{format_substituent_prefixes(grouped)}{parent}({count + 2})"


def lewis_adduct_mol(mol):
    """The mol split into its neutral donor and acceptor components, or None."""
    editable = Chem.RWMol(mol)
    cuts = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() == Chem.BondType.DATIVE:
            if a.GetAtomicNum() in _DONORS and b.GetAtomicNum() in _ACCEPTORS:
                cuts.append((a.GetIdx(), b.GetIdx()))
            continue
        for donor, acceptor in ((a, b), (b, a)):
            if (
                donor.GetFormalCharge() == 1
                and acceptor.GetFormalCharge() == -1
                and donor.GetAtomicNum() in _DONORS
                and acceptor.GetAtomicNum() in _ACCEPTORS
                and bond.GetBondTypeAsDouble() == 1.0
            ):
                cuts.append((donor.GetIdx(), acceptor.GetIdx()))
    if not cuts or any(a.GetFormalCharge() and not any(a.GetIdx() in cut for cut in cuts) for a in mol.GetAtoms()):
        return None
    if any(b.GetBondType() == Chem.BondType.DATIVE for b in mol.GetBonds()) and sum(
        b.GetBondType() == Chem.BondType.DATIVE for b in mol.GetBonds()
    ) != len(cuts):
        return None
    if any(a.GetAtomicNum() > 20 and a.GetAtomicNum() not in _ACCEPTORS and a.GetAtomicNum() not in (33, 34, 35, 53) for a in mol.GetAtoms()):
        return None
    for donor, acceptor in cuts:
        editable.RemoveBond(donor, acceptor)
        if len(cuts) == 1:
            editable.GetAtomWithIdx(donor).SetBoolProp("_adduct_donor", True)
            editable.GetAtomWithIdx(acceptor).SetBoolProp("_adduct_acceptor", True)
        for index in (donor, acceptor):
            atom = editable.GetAtomWithIdx(index)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
    parts = editable.GetMol()
    try:
        Chem.SanitizeMol(parts)
    except Exception:
        return None
    return parts if len(Chem.GetMolFrags(parts)) > 1 else None
