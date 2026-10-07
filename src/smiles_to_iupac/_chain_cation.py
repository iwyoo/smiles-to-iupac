"""Cationic centres of an acyclic chain of oxygen and sulfur named by skeletal replacement (P-73.4, P-51.4): the
neutral skeleton takes the 'oxa'/'thia' prefixes and the cationic centres follow as 'ium' with the lowest locants after
the heteroatoms ('3,12-dimethyl-6,9-dioxa-3,12-dithiatetradecane-3,12-diium'); substituents on the chain are prefixes.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents, substituent_locant_set_and_citation
from ._numerals import alkane_name, multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_A_PREFIX = {8: "oxa", 16: "thia"}
_SENIORITY = [8, 16]
_SIGMA = {8: 3, 16: 3}
_MIN_HETEROUNITS = 4


def has_chain_cation_shape(mol) -> bool:
    centres = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if not centres or len(Chem.GetMolFrags(mol)) != 1 or mol.GetRingInfo().NumRings():
        return False
    if any(a.GetAtomicNum() not in (6, 8, 16) or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return False
    if any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()):
        return False
    if any(a.GetAtomicNum() == 6 and a.GetFormalCharge() for a in mol.GetAtoms()):
        return False
    if any(a.GetFormalCharge() != 1 or a.GetDegree() + a.GetTotalNumHs() != _SIGMA[a.GetAtomicNum()] for a in centres):
        return False
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 6) >= _MIN_HETEROUNITS


def _paths(graph, leaves):
    found = []
    for start in leaves:
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            onward = [n for n in graph[node] if n not in path]
            if not onward and len(path) > 1:
                found.append(path)
            for n in onward:
                stack.append((n, path + [n]))
    return found


def name_chain_cation(mol) -> str:
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        raise UnsupportedStructure("stereodescriptors of a skeletal replacement cation are not supported yet")
    graph = adjacency(mol)
    hetero = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 6}
    leaves = [a for a, n in graph.items() if len(n) == 1]
    candidates = [
        p
        for p in _paths(graph, leaves)
        if hetero <= set(p)
        and mol.GetAtomWithIdx(p[0]).GetAtomicNum() == 6
        and mol.GetAtomWithIdx(p[-1]).GetAtomicNum() == 6
        and p[0] < p[-1]
    ]
    if not candidates:
        raise UnsupportedStructure("the heteroatoms do not lie on one chain ending in carbon")
    longest = max(len(p) for p in candidates)
    chain = next(p for p in candidates if len(p) == longest)
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and mol.GetAtomWithIdx(b).GetAtomicNum() != 6 for a, b in zip(chain, chain[1:])):
        raise UnsupportedStructure("adjacent heteroatoms in a skeletal replacement chain are not supported here")
    halogens = halogen_substituents(mol)
    chain_set = set(chain)
    best = None
    for ordered in (chain, chain[::-1]):
        position = {atom: i + 1 for i, atom in enumerate(ordered)}
        hetero_locants = {z: [position[a] for a in ordered if mol.GetAtomWithIdx(a).GetAtomicNum() == z] for z in _SENIORITY}
        ium = [position[a] for a in ordered if mol.GetAtomWithIdx(a).GetFormalCharge()]
        entries = {}
        for atom in ordered:
            for n in graph[atom]:
                if n not in chain_set:
                    entries.setdefault(position[atom], []).append(name_branch(graph, n, atom, halogens, mol=mol))
        grouped = group_substituents(entries)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        key = (
            sorted(l for v in hetero_locants.values() for l in v),
            [hetero_locants[z] for z in _SENIORITY],
            ium,
            locant_set,
            citation,
        )
        if best is None or key < best[0]:
            best = (key, hetero_locants, ium, grouped)
    _, hetero_locants, ium, grouped = best
    parts = [
        f"{','.join(map(str, hetero_locants[z]))}-{multiplying_prefix(len(hetero_locants[z])) if len(hetero_locants[z]) > 1 else ''}{_A_PREFIX[z]}"
        for z in _SENIORITY
        if hetero_locants[z]
    ]
    parent = alkane_name(len(chain))
    if len(ium) == 1:
        ending = f"{parent[:-1]}-{ium[0]}-ium"
    else:
        ending = f"{parent}-{','.join(map(str, ium))}-{multiplying_prefix(len(ium))}ium"
    core = "-".join(parts) + ending
    prefixes = format_substituent_prefixes(grouped) if grouped else ""
    return f"{prefixes}-{core}" if prefixes else core
