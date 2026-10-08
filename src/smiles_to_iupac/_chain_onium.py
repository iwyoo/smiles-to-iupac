"""Onium cations of unbranched homonuclear heteroatom chains (P-73.1.1, P-73.1.2): hydrazin-1-ium, hydrazine-1,2-diium,
1-methylhydrazin-1-ium, 2-benzoyldioxidan-1-ium. The cationic atoms of the chain take the lowest locants, then the
substituents."""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    substituent_locant_set_and_citation,
)
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_STEM = {7: ("azane", 3), 8: ("oxidane", 2), 15: ("phosphane", 3), 16: ("sulfane", 2), 34: ("selane", 2), 52: ("tellane", 2)}


def _chain(mol, start):
    """The unbranched chain of single-bonded atoms of the element of `start`, as an ordered list, or None."""
    z = mol.GetAtomWithIdx(start).GetAtomicNum()
    component = {start}
    stack = [start]
    while stack:
        current = stack.pop()
        for n in mol.GetAtomWithIdx(current).GetNeighbors():
            if (
                n.GetAtomicNum() == z
                and n.GetIdx() not in component
                and mol.GetBondBetweenAtoms(current, n.GetIdx()).GetBondTypeAsDouble() == 1.0
            ):
                component.add(n.GetIdx())
                stack.append(n.GetIdx())
    if any(mol.GetAtomWithIdx(a).IsInRing() for a in component):
        return None
    ends = [a for a in component if sum(n.GetIdx() in component for n in mol.GetAtomWithIdx(a).GetNeighbors()) <= 1]
    if len(component) < 2 or len(ends) != 2 or any(
        sum(n.GetIdx() in component for n in mol.GetAtomWithIdx(a).GetNeighbors()) > 2 for a in component
    ):
        return None
    order = [ends[0]]
    while len(order) < len(component):
        order.append(
            next(n.GetIdx() for n in mol.GetAtomWithIdx(order[-1]).GetNeighbors() if n.GetIdx() in component and n.GetIdx() not in order)
        )
    return order


def _match(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if not charged or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() != 1 for a in charged):
        return None
    z = charged[0].GetAtomicNum()
    if z not in _STEM or any(a.GetAtomicNum() != z for a in charged) or any(a.GetIsotope() for a in mol.GetAtoms()):
        return None
    chain = _chain(mol, charged[0].GetIdx())
    if chain is None or not {a.GetIdx() for a in charged} <= set(chain):
        return None
    standard = _STEM[z][1]
    for atom in charged:
        if sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs() != standard + 1:
            return None
    if any(a.GetFormalCharge() == 0 and a.GetAtomicNum() == z and a.GetIdx() in chain and a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    return z, chain, {a.GetIdx() for a in charged}


def has_chain_onium_shape(mol) -> bool:
    return _match(mol) is not None


def name_chain_onium(mol) -> str:
    found = _match(mol)
    if found is None:
        raise UnsupportedStructure("not an onium cation of a heteroatom chain")
    z, chain, cations = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    hydrazide = _hydrazide_ium(mol, graph, z, chain, cations, halogens, aromatic)
    if hydrazide is not None:
        return hydrazide
    chain_set = set(chain)
    best = None
    for walk in (chain, chain[::-1]):
        positions = {a: i + 1 for i, a in enumerate(walk)}
        subs = {}
        for atom in walk:
            for n in graph[atom]:
                if n in chain_set:
                    continue
                subs.setdefault(positions[atom], []).append(
                    name_branch(graph, n, atom, halogens, aromatic, mol=mol, unsaturated=True)
                )
        grouped = group_substituents(subs)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        key = (sorted(positions[a] for a in cations), locant_set, citation)
        if best is None or key < best[0]:
            best = (key, grouped)
    key, grouped = best
    count = len(cations)
    n = len(chain)
    stem = "hydrazin" if z == 7 and n == 2 else multiplying_prefix(n) + _STEM[z][0][:-1]
    locants = ",".join(map(str, key[0]))
    if count == 1:
        parent = f"{stem}-{locants}-ium"
    else:
        parent = f"{stem}e-{locants}-{multiplying_prefix(count)}ium"
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    return prefix + parent


def _hydrazide_ium(mol, graph, z, chain, cations, halogens, aromatic):
    """P-73.1.2.1: an acylated hydrazine whose other nitrogen is the cationic centre is a hydrazid-N'-ium of the acid,
    N',N',N'-trimethylbenzohydrazid-N'-ium, not a hydrazin-1-ium with an acyl prefix."""
    from .core import smiles_to_iupac

    if z != 7 or len(chain) != 2 or len(cations) != 1:
        return None
    (cation,) = cations
    acylated = next(a for a in chain if a != cation)
    carbonyl = [
        n for n in graph[acylated]
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 6
        and any(
            m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n, m.GetIdx()).GetBondTypeAsDouble() == 2.0
            for m in mol.GetAtomWithIdx(n).GetNeighbors()
        )
    ]
    if len(carbonyl) != 1:
        return None
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        if atom.GetIdx() in chain:
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(2 if atom.GetIdx() == cation else 1)
    branches = {a: [n for n in graph[a] if n not in chain and n != carbonyl[0]] for a in chain}
    removed = set()
    for a in chain:
        for n in branches[a]:
            stack, seen = [n], set()
            while stack:
                current = stack.pop()
                if current in seen:
                    continue
                seen.add(current)
                stack.extend(m for m in graph[current] if m not in chain and m != carbonyl[0] and m not in seen)
            removed |= seen
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    neutral = editable.GetMol()
    try:
        Chem.SanitizeMol(neutral)
        parent = smiles_to_iupac(Chem.MolToSmiles(neutral))
    except UnsupportedStructure:
        return None
    if not parent.endswith("hydrazide") or " " in parent or "-" in parent[: parent.rfind("hydrazide")] or parent.startswith("("):
        return None
    if parent.count("(") or any(ch.isdigit() for ch in parent):
        return None
    entries = {}
    for atom, label in ((acylated, "N"), (cation, "N'")):
        for n in branches[atom]:
            entries.setdefault(label, []).append(name_branch(graph, n, atom, halogens, aromatic, mol=mol, unsaturated=True))
    grouped = group_substituents(entries)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    return prefix + parent[:-1] + "-N'-ium"
