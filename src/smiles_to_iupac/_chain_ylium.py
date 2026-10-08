"""Ylium cations of unbranched homonuclear heteroatom chains (P-73.2.2.1.2, P-73.2.3.3, P-73.2.3.4): hydrazin-1-ylium,
hydrazine-1,1-bis(ylium), triaz-1-en-1-ylium, heptamethyltrisilan-2-ylium, tert-butylperoxylium, phenyldisulfanylium.
The cationic atoms take the lowest locants, then the double bonds, then the substituents."""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    substituent_locant_set_and_citation,
)
from ._locant_omission import omits_all_locants
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_STEM = {
    7: ("az", 3), 8: ("oxid", 2), 14: ("sil", 4), 15: ("phosph", 3), 16: ("sulf", 2), 33: ("ars", 3), 34: ("sel", 2),
    51: ("stib", 3), 52: ("tellur", 2),
}
_GROUP = {8: "perox", 16: "disulfan", 34: "diselan", 52: "ditellan"}
_BIS = {2: "bis", 3: "tris", 4: "tetrakis"}


def _chain(mol, start):
    z = mol.GetAtomWithIdx(start).GetAtomicNum()
    component = {start}
    stack = [start]
    while stack:
        current = stack.pop()
        for n in mol.GetAtomWithIdx(current).GetNeighbors():
            bond = mol.GetBondBetweenAtoms(current, n.GetIdx())
            if n.GetAtomicNum() == z and n.GetIdx() not in component and bond.GetBondTypeAsDouble() in (1.0, 2.0):
                component.add(n.GetIdx())
                stack.append(n.GetIdx())
    degree = lambda a: sum(n.GetIdx() in component for n in mol.GetAtomWithIdx(a).GetNeighbors())
    ends = [a for a in component if degree(a) <= 1]
    if len(component) < 2 or len(ends) != 2 or any(degree(a) > 2 for a in component):
        return None
    if any(mol.GetAtomWithIdx(a).IsInRing() for a in component):
        return None
    order = [ends[0]]
    while len(order) < len(component):
        order.append(
            next(n.GetIdx() for n in mol.GetAtomWithIdx(order[-1]).GetNeighbors() if n.GetIdx() in component and n.GetIdx() not in order)
        )
    return order


def _match(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if not charged or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() not in (1, 2) for a in charged):
        return None
    z = charged[0].GetAtomicNum()
    if z not in _STEM or any(a.GetAtomicNum() != z for a in charged) or any(a.GetIsotope() for a in mol.GetAtoms()):
        return None
    for atom in charged:
        bonds = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()
        if bonds != _STEM[z][1] - atom.GetFormalCharge():
            return None
    chain = _chain(mol, charged[0].GetIdx())
    if chain is None or not {a.GetIdx() for a in charged} <= set(chain):
        return None
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms() if a.GetFormalCharge() == 0):
        return None
    return z, chain, {a.GetIdx(): a.GetFormalCharge() for a in charged}


def has_chain_ylium_shape(mol) -> bool:
    return _match(mol) is not None


def name_chain_ylium(mol) -> str:
    found = _match(mol)
    if found is None:
        raise UnsupportedStructure("not a ylium cation of a heteroatom chain")
    z, chain, cations = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    chain_set = set(chain)
    best = None
    for walk in (chain, chain[::-1]):
        positions = {a: i + 1 for i, a in enumerate(walk)}
        enes = [
            i + 1
            for i in range(len(walk) - 1)
            if mol.GetBondBetweenAtoms(walk[i], walk[i + 1]).GetBondTypeAsDouble() == 2.0
        ]
        subs = {}
        for atom in walk:
            for n in graph[atom]:
                if n not in chain_set:
                    subs.setdefault(positions[atom], []).append(
                        name_branch(graph, n, atom, halogens, aromatic, mol=mol, unsaturated=True)
                    )
        grouped = group_substituents(subs)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        ylium = sorted(p for a, c in cations.items() for p in [positions[a]] * c)
        key = (ylium, enes, locant_set, citation)
        if best is None or key < best[0]:
            best = (key, grouped)
    (ylium, enes, _, _), grouped = best
    n = len(chain)
    omitted = z != 7 and omits_all_locants(mol, set(chain), grouped, single_kind=False)
    prefix = format_substituent_prefixes(grouped, omit_all=omitted) if grouped else ""
    root = _STEM[z][0]
    if n == 2 and z != 7 and not enes and len(ylium) == 1:
        group = _GROUP[z] if grouped or z != 8 else "dioxidan"
        return prefix.removeprefix("2-") + group + "ylium"
    if n == 2 and enes and len(ylium) == 1:
        return prefix + multiplying_prefix(n) + root + "enylium"
    if enes:
        ene = "en" if len(enes) == 1 else multiplying_prefix(len(enes)) + "en"
        stem = f"{multiplying_prefix(n)}{root}-{','.join(map(str, enes))}-{ene}"
    else:
        stem = "hydrazin" if z == 7 and n == 2 else f"{multiplying_prefix(n)}{root}an"
    locants = ",".join(map(str, ylium))
    if len(ylium) == 1:
        return f"{prefix}{stem}-{locants}-ylium"
    return f"{prefix}{stem}e-{locants}-{_BIS[len(ylium)]}(ylium)"
