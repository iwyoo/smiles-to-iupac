"""Acyclic polyboranes with boron-boron bonds, and Lewis adducts of boranes (P-68.1.1.2.1, P-68.1.6).

An unbranched chain of n boron atoms is the multiplied parent hydride 'diborane', 'triborane', ...; the number of
hydrogen atoms of the parent hydride, n + 2 for the chain of three-coordinate borons, follows in parentheses. A donor
bonded to an acceptor by a dative bond, or drawn as the charge-separated N(+)-B(-) pair, is named as the adduct of the
two neutral components (P-68.1.6.2), which is not a preferred IUPAC name.
"""

from rdkit import Chem

from ._cited_group import cited_group
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, numerical_term
from ._substituents import format_substituent_prefixes

_ACCEPTORS = {5, 13, 31}
_DONORS = {7, 8, 15, 16}


def _boron_path(mol, borons):
    graph = {b: [n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetAtomicNum() == 5] for b in borons}
    if any(len(v) > 2 for v in graph.values()):
        return None
    ends = [b for b, v in graph.items() if len(v) == 1]
    if len(ends) != 2:
        return None
    path = [ends[0]]
    while len(path) < len(borons):
        following = [n for n in graph[path[-1]] if n not in path]
        if len(following) != 1:
            return None
        path.append(following[0])
    return path


def polyborane_name(mol):
    borons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 5]
    if len(borons) < 2 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(mol.GetAtomWithIdx(b).GetFormalCharge() or mol.GetAtomWithIdx(b).GetIsotope() for b in borons):
        return None
    if any(bond.GetBondTypeAsDouble() != 1.0 for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 5 and bond.GetEndAtom().GetAtomicNum() == 5):
        return None
    path = _boron_path(mol, borons)
    if path is None or any(mol.GetAtomWithIdx(b).GetDegree() + mol.GetAtomWithIdx(b).GetTotalNumHs() != 3 for b in path):
        return None
    graph = adjacency(mol)
    entries = []
    try:
        for position, boron in enumerate(path):
            for n in graph[boron]:
                if mol.GetAtomWithIdx(n).GetAtomicNum() != 5:
                    entries.append((position, cited_group(mol, graph, n, boron)))
    except UnsupportedStructure:
        return None
    count = len(path)
    best = None
    for reverse in (False, True):
        locants = [(count - position if reverse else position + 1, name, compound) for position, (name, compound) in entries]
        key = (sorted(l for l, _, _ in locants), [l for l, _, _ in sorted(locants, key=lambda item: alpha_sort_key(item[1]))])
        if best is None or key < best[0]:
            best = (key, locants)
    grouped = {}
    for locant, name, compound in best[1]:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    for info in grouped.values():
        info["locants"].sort()
    return f"{format_substituent_prefixes(grouped)}{numerical_term(count)}borane({count + 2})"


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
