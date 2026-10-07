"""Cationic centres of an acyclic chain of chalcogen and pnictogen heterounits named by skeletal replacement (P-73.4,
P-51.4): the neutral skeleton takes the 'oxa'/'thia'/'phospha' prefixes, a heterounit with a nonstandard bonding number
carries its lambda number, and the cationic centres follow as 'ium' with the lowest locants after the heteroatoms
('3,12-dimethyl-6,9-dioxa-3,12-dithiatetradecane-3,12-diium'); substituents on the chain are prefixes.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents, substituent_locant_set_and_citation
from ._numerals import alkane_name, multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_A_PREFIX = {8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 15: "phospha", 33: "arsa"}
_SENIORITY = [8, 16, 34, 52, 15, 33]
_SIGMA = {8: 3, 16: 3, 34: 3, 52: 3, 15: 4, 33: 4}
_STANDARD = {8: 2, 16: 2, 34: 2, 52: 2, 15: 3, 33: 3}
_OXO = {8, 16, 34, 52}
_MIN_HETEROUNITS = 4


def _terminal_chalcogens(mol):
    """Doubly bonded terminal chalcogens: the oxo-type substituents of a heterounit with a nonstandard bonding number."""
    found = set()
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        for end, other in ((bond.GetBeginAtom(), bond.GetEndAtom()), (bond.GetEndAtom(), bond.GetBeginAtom())):
            if end.GetAtomicNum() in _OXO and end.GetDegree() == 1 and not end.GetFormalCharge() and other.GetAtomicNum() in _STANDARD:
                found.add(end.GetIdx())
    return found


def _chain_heteroatoms(mol):
    return {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 6} - _terminal_chalcogens(mol)


def _valence(atom):
    return int(round(atom.GetTotalValence()))


def has_chain_cation_shape(mol) -> bool:
    centres = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if not centres or len(Chem.GetMolFrags(mol)) != 1:
        return False
    if any(a.IsInRing() and (a.GetAtomicNum() != 6 or a.GetFormalCharge()) for a in mol.GetAtoms()):
        return False
    if any(a.GetAtomicNum() not in (6, *_STANDARD) or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return False
    terminal = _terminal_chalcogens(mol)
    if any(
        b.GetBondTypeAsDouble() != 1.0 and not b.IsInRing() and not (b.GetBeginAtomIdx() in terminal or b.GetEndAtomIdx() in terminal)
        for b in mol.GetBonds()
    ):
        return False
    if any(a.GetAtomicNum() == 6 and a.GetFormalCharge() for a in mol.GetAtoms()):
        return False
    if any(a.GetFormalCharge() != 1 or a.GetAtomicNum() not in _SIGMA or _valence(a) != _SIGMA[a.GetAtomicNum()] for a in centres):
        return False
    chain = _chain_heteroatoms(mol)
    if any(
        not a.GetFormalCharge() and (_valence(a) < _STANDARD[a.GetAtomicNum()] or (_valence(a) - _STANDARD[a.GetAtomicNum()]) % 2)
        for a in (mol.GetAtomWithIdx(i) for i in chain)
    ):
        return False
    return len(chain) >= _MIN_HETEROUNITS


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
    hetero = _chain_heteroatoms(mol)
    terminal = _terminal_chalcogens(mol)
    ring_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
    skeleton = {a: [n for n in ns if n not in ring_atoms and n not in terminal] for a, ns in graph.items() if a not in ring_atoms and a not in terminal}
    leaves = [a for a, n in skeleton.items() if len(n) == 1]
    candidates = [
        p
        for p in _paths(skeleton, leaves)
        if hetero <= set(p)
        and mol.GetAtomWithIdx(p[0]).GetAtomicNum() == 6
        and mol.GetAtomWithIdx(p[-1]).GetAtomicNum() == 6
        and p[0] < p[-1]
    ]
    if not candidates:
        raise UnsupportedStructure("the heteroatoms do not lie on one chain ending in carbon")
    longest = max(len(p) for p in candidates)
    chain = next(p for p in candidates if len(p) == longest)
    halogens = halogen_substituents(mol)
    chain_set = set(chain)
    nonstandard = {
        a: _valence(mol.GetAtomWithIdx(a)) - mol.GetAtomWithIdx(a).GetFormalCharge()
        for a in chain
        if a in hetero and _valence(mol.GetAtomWithIdx(a)) - mol.GetAtomWithIdx(a).GetFormalCharge() > _STANDARD[mol.GetAtomWithIdx(a).GetAtomicNum()]
    }
    best = None
    for ordered in (chain, chain[::-1]):
        position = {atom: i + 1 for i, atom in enumerate(ordered)}
        hetero_locants = {z: [position[a] for a in ordered if mol.GetAtomWithIdx(a).GetAtomicNum() == z] for z in _SENIORITY}
        ium = [position[a] for a in ordered if mol.GetAtomWithIdx(a).GetFormalCharge()]
        lam = {position[a]: n for a, n in nonstandard.items()}
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
            sorted((-n, l) for l, n in lam.items()),
            ium,
            locant_set,
            citation,
        )
        if best is None or key < best[0]:
            best = (key, hetero_locants, ium, grouped, lam)
    _, hetero_locants, ium, grouped, lam = best
    parts = [
        f"{','.join(f'{l}λ{lam[l]}' if l in lam else str(l) for l in hetero_locants[z])}-{multiplying_prefix(len(hetero_locants[z])) if len(hetero_locants[z]) > 1 else ''}{_A_PREFIX[z]}"
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
