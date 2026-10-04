"""Unbranched chains of one heteroatom (disilane, triazane, 1,2-dimethyldiphosphane) as parent hydrides
(P-21.2.1, P-44.1.2.2): Group 14 and 15 atoms outrank carbon, so the chain is the parent and carbon groups are prefixes."""

from ._common import (
    UnsupportedStructure,
    group_substituents,
    substituent_locant_set_and_citation,
)
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_STEMS = {7: "azane", 14: "silane", 15: "phosphane", 32: "germane", 33: "arsane", 50: "stannane", 51: "stibane", 82: "plumbane", 83: "bismuthane"}


def _chain_atoms(mol, graph):
    elements = {a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS and not a.IsInRing()}
    if len(elements) != 1:
        return None
    (z,) = elements
    atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == z}
    if any(mol.GetAtomWithIdx(a).IsInRing() or mol.GetAtomWithIdx(a).GetFormalCharge() for a in atoms):
        return None
    ends = [a for a in atoms if sum(n in atoms for n in graph[a]) <= 1]
    if len(atoms) < (3 if z == 7 else 2) or len(ends) != 2 or any(sum(n in atoms for n in graph[a]) > 2 for a in atoms):
        return None
    start = ends[0]
    chain, previous = [start], None
    while True:
        nxt = [n for n in graph[chain[-1]] if n in atoms and n != previous]
        if not nxt:
            break
        previous = chain[-1]
        chain.append(nxt[0])
    if len(chain) != len(atoms):
        return None
    for a, b in zip(chain, chain[1:]):
        if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0:
            return None
    return z, chain


def _siloxane_chain(mol, graph):
    atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (8, 14) and not a.IsInRing()}
    silicon = [a for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 14]
    oxygen = [a for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]
    if len(silicon) < 2 or len(oxygen) != len(silicon) - 1:
        return None
    if any(sum(n in atoms for n in graph[a]) > 2 for a in atoms):
        return None
    if any(sum(n in atoms for n in graph[a]) != 2 for a in oxygen) or any(mol.GetAtomWithIdx(a).GetDegree() != 2 for a in oxygen):
        return None
    ends = [a for a in silicon if sum(n in atoms for n in graph[a]) == 1]
    if len(ends) != 2:
        return None
    chain, previous = [ends[0]], None
    while True:
        nxt = [n for n in graph[chain[-1]] if n in atoms and n != previous]
        if not nxt:
            break
        previous = chain[-1]
        chain.append(nxt[0])
    if len(chain) != len(atoms) or any(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0 for a, b in zip(chain, chain[1:])):
        return None
    return chain


def name_hydride_chain(mol, graph, halogens, aromatic_atoms):
    found = _chain_atoms(mol, graph)
    siloxane = None
    if found is None:
        siloxane = _siloxane_chain(mol, graph)
        if siloxane is None:
            return None
        z, chain = 14, siloxane
    else:
        z, chain = found
    chain_set = set(chain)
    best = None
    for candidate in (chain, chain[::-1]):
        substituents = {}
        for i, atom in enumerate(candidate):
            for n in graph[atom]:
                if n in chain_set:
                    continue
                if mol.GetBondBetweenAtoms(atom, n).GetBondTypeAsDouble() != 1.0:
                    raise UnsupportedStructure("a multiple bond on a heteroatom chain is not supported yet")
                name, compound = name_branch(graph, n, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                substituents.setdefault(i + 1, []).append((name, compound))
        grouped = group_substituents(substituents)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        omit = len(chain) == 2 and sum(len(info["locants"]) for info in grouped.values()) == 1
        if siloxane is not None:
            parent = f"{multiplying_prefix((len(chain) + 1) // 2)}siloxane"
            omit = False
        else:
            parent = f"{multiplying_prefix(len(chain))}{_STEMS[z]}"
        name = format_substituent_prefixes(grouped, omit_locants=omit) + parent
        key = (locant_set, citation, name)
        if best is None or key < best[0]:
            best = (key, name)
    return best[1]
