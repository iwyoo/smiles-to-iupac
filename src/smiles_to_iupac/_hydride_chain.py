"""Unbranched chains of one heteroatom (disilane, triazane, 1,2-dimethyldiphosphane) as parent hydrides
(P-21.2.1, P-44.1.2.2): Group 14 and 15 atoms outrank carbon, so the chain is the parent and carbon groups are prefixes."""

from ._common import (
    UnsupportedStructure,
    group_substituents,
    nonstandard_bonding,
    substituent_locant_set_and_citation,
)
from ._numerals import multiplying_prefix
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_STEMS = {7: "azane", 14: "silane", 15: "phosphane", 32: "germane", 33: "arsane", 50: "stannane", 51: "stibane", 82: "plumbane", 83: "bismuthane"}


def _is_nitrile_nitrogen(atom):
    return (
        atom.GetAtomicNum() == 7
        and atom.GetDegree() == 1
        and atom.GetBonds()[0].GetBondTypeAsDouble() == 3.0
        and atom.GetNeighbors()[0].GetAtomicNum() == 6
    )


def _chain_atoms(mol, graph, skip_nitrogen=False, allow_double=False):
    elements = {
        a.GetAtomicNum()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in _STEMS
        and not a.IsInRing()
        and not (skip_nitrogen and a.GetAtomicNum() == 7)
        and not (allow_double and _is_nitrile_nitrogen(a))
    }
    if len(elements) != 1:
        return None
    (z,) = elements
    atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == z and not (allow_double and _is_nitrile_nitrogen(a))}
    if any(mol.GetAtomWithIdx(a).IsInRing() or mol.GetAtomWithIdx(a).GetFormalCharge() for a in atoms):
        return None
    ends = [a for a in atoms if sum(n in atoms for n in graph[a]) <= 1]
    if len(atoms) < (3 if z == 7 and not allow_double else 2) or len(ends) != 2 or any(sum(n in atoms for n in graph[a]) > 2 for a in atoms):
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
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order != 1.0 and not (allow_double and order == 2.0):
            return None
    if z == 7 and len(chain) == 2 and mol.GetBondBetweenAtoms(chain[0], chain[1]).GetBondTypeAsDouble() != 2.0:
        return None
    return z, chain


# P-21.2.3.1: seniority O > S > Se > Te > P > ... > Tl; N is excluded because amine names are preferred
_ALTERNATING_ORDER = (8, 16, 34, 52, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81)
_A_TERMS = {
    8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma", 14: "sila",
    32: "germa", 50: "stanna", 82: "plumba", 5: "bora", 13: "alumina", 31: "gallia", 49: "india", 81: "thallia",
}
_DIVALENT = {8, 16, 34, 52}


def _alternating_chain(mol, graph):
    """(terminal element, inner element, chain) of an unbranched a(ba)n chain of two heteroatoms (P-21.2.3.1)."""
    atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in _ALTERNATING_ORDER and not a.IsInRing()}
    elements = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in atoms}
    if len(elements) != 2:
        return None
    inner_z, terminal_z = sorted(elements, key=_ALTERNATING_ORDER.index)
    terminals = [a for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == terminal_z]
    inner = [a for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == inner_z]
    if len(terminals) < 2 or len(inner) != len(terminals) - 1:
        return None
    if any(sum(n in atoms for n in graph[a]) > 2 for a in atoms):
        return None
    if any(sum(n in atoms for n in graph[a]) != 2 for a in inner):
        return None
    if inner_z in _DIVALENT and any(mol.GetAtomWithIdx(a).GetDegree() != 2 for a in inner):
        return None
    ends = [a for a in terminals if sum(n in atoms for n in graph[a]) == 1]
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
    return terminal_z, inner_z, chain


def _alternating_parent(terminal_z, inner_z, terminal_count):
    terminal, inner = _A_TERMS[terminal_z], _A_TERMS[inner_z]
    if inner[0] in "aeiou":
        terminal = terminal[:-1]
    return multiplying_prefix(terminal_count) + terminal + inner[:-1] + "ane"


def name_hydride_chain(mol, graph, halogens, aromatic_atoms):
    found = _chain_atoms(mol, graph)
    alternating = None
    if found is None:
        alternating = _alternating_chain(mol, graph)
        if alternating is None:
            return None
        chain = alternating[2]
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
                if mol.GetBondBetweenAtoms(atom, n).GetBondTypeAsDouble() not in (1.0, 2.0, 3.0) or (
                    mol.GetBondBetweenAtoms(atom, n).GetBondTypeAsDouble() != 1.0 and mol.GetAtomWithIdx(n).GetAtomicNum() != 6
                ):
                    raise UnsupportedStructure("a multiple bond to a non-carbon group on a heteroatom chain is not supported yet")
                name, compound = name_branch(graph, n, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                substituents.setdefault(i + 1, []).append((name, compound))
        grouped = group_substituents(substituents)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        lam = {i + 1: n for i, atom in enumerate(candidate) if (n := nonstandard_bonding(mol.GetAtomWithIdx(atom)))}
        omit = not lam and (
            (len(chain) == 2 and sum(len(info["locants"]) for info in grouped.values()) == 1)
            or all(mol.GetAtomWithIdx(a).GetTotalNumHs() == 0 for a in chain)
        )
        if alternating is not None:
            parent = _alternating_parent(alternating[0], alternating[1], (len(chain) + 1) // 2)
            omit = not lam and len(chain) == 3 and sum(len(info["locants"]) for info in grouped.values()) == 1
        else:
            parent = f"{multiplying_prefix(len(chain))}{_STEMS[z]}"
        lam_text = ",".join(f"{p}\u03bb{lam[p]}" for p in sorted(lam))
        if lam_text:
            parent = f"{lam_text}-{parent}"
        if omit and len(grouped) > 1:
            entries = [(n, info["compound"]) for n, info in grouped.items() for _ in info["locants"]]
            name = format_mononuclear_prefixes(entries) + parent
        else:
            prefixes = format_substituent_prefixes(grouped, omit_locants=omit)
            name = prefixes + ("-" if prefixes and lam_text else "") + parent
        key = (sorted(lam), [-lam[p] for p in sorted(lam)], locant_set, citation, name)
        if best is None or key < best[0]:
            best = (key, name)
    return best[1]
