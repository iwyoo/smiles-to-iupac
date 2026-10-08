"""Unbranched chains of one heteroatom (disilane, triazane, 1,2-dimethyldiphosphane) as parent hydrides
(P-21.2.1, P-44.1.2.2): Group 14 and 15 atoms outrank carbon, so the chain is the parent and carbon groups are prefixes."""

from ._common import (
    UnsupportedStructure,
    unsaturation_suffix,
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


def _is_nitro_or_nitroso(atom):
    """A nitro or nitroso nitrogen: bonded to one non-oxygen atom and otherwise only oxygen (a substituent, not a chain atom)."""
    if atom.GetAtomicNum() != 7:
        return False
    others = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 8]
    mol = atom.GetOwningMol()
    return len(others) == 1 and any(
        n.GetAtomicNum() == 8
        and (mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0 or n.GetFormalCharge() == -1)
        for n in atom.GetNeighbors()
    )


def _is_pseudohalide_nitrogen(atom):
    """The N of an isocyanato-type group N=C=X (X = O, S, Se, Te): a compulsory prefix, not a chain atom (P-58.3.2)."""
    if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetTotalNumHs() or atom.GetFormalCharge() or atom.IsInRing():
        return False
    mol = atom.GetOwningMol()
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    return len(carbons) == 1 and carbons[0].GetDegree() == 2 and any(
        n.GetAtomicNum() in (8, 16, 34, 52) and n.GetIdx() != atom.GetIdx()
        and mol.GetBondBetweenAtoms(carbons[0].GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in carbons[0].GetNeighbors()
    )


def _chain_atoms(mol, graph, skip_nitrogen=False, allow_double=False, allow_triple=False):
    elements = {
        a.GetAtomicNum()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in _STEMS
        and not _is_nitro_or_nitroso(a)
        and not _is_pseudohalide_nitrogen(a)
        and not a.IsInRing()
        and not (skip_nitrogen and a.GetAtomicNum() == 7)
        and not (allow_double and _is_nitrile_nitrogen(a))
    }
    if len(elements) != 1:
        return None
    (z,) = elements
    atoms = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == z
        and not (allow_double and _is_nitrile_nitrogen(a))
        and not _is_nitro_or_nitroso(a)
        and not _is_pseudohalide_nitrogen(a)
    }
    if any(mol.GetAtomWithIdx(a).IsInRing() or mol.GetAtomWithIdx(a).GetFormalCharge() for a in atoms):
        return None
    ends = [a for a in atoms if sum(n in atoms for n in graph[a]) <= 1]
    pseudohalide = any(_is_pseudohalide_nitrogen(a) for a in mol.GetAtoms())
    minimum = 3 if z == 7 and not allow_double and not pseudohalide else 2
    if len(atoms) < minimum or len(ends) != 2 or any(sum(n in atoms for n in graph[a]) > 2 for a in atoms):
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
        if order != 1.0 and not (allow_double and order == 2.0) and not (allow_triple and order == 3.0):
            return None
    if z == 7 and len(chain) == 2 and mol.GetBondBetweenAtoms(chain[0], chain[1]).GetBondTypeAsDouble() != 2.0 and not pseudohalide:
        return None
    return z, chain


# P-21.2.3.1: seniority O > S > Se > Te > P > ... > Tl; N is excluded because amine names are preferred
_ALTERNATING_ORDER = (8, 16, 34, 52, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81)
_A_TERMS = {
    8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma", 14: "sila",
    32: "germa", 50: "stanna", 82: "plumba", 5: "bora", 13: "alumina", 31: "galla", 49: "inda", 81: "thalla",
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
    if len(chain) != len(atoms):
        return None
    return terminal_z, inner_z, chain


def _alternating_parent(terminal_z, inner_z, terminal_count, ene=(), yne=()):
    terminal, inner = _A_TERMS[terminal_z], _A_TERMS[inner_z]
    if inner[0] in "aeiou":
        terminal = terminal[:-1]
    stem = multiplying_prefix(terminal_count) + terminal + inner[:-1]
    if not ene and not yne:
        return stem + "ane"
    body, needs_a = unsaturation_suffix(ene, yne)
    return f"{stem}{'a' if needs_a else ''}-{body}"


def _split_chain(mol, graph):
    """(element, chain) of the senior unbranched run of one Group 14/15 element when other atoms (a disulfanyl link)
    join several such runs: the longest run, then the one carrying more substituents (P-44.1.2.2, P-44.3)."""
    elements = {
        a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS and not a.IsInRing() and not _is_nitro_or_nitroso(a)
    }
    if len(elements) != 1 or any(a.GetFormalCharge() for a in mol.GetAtoms()):
        return None
    (z,) = elements
    atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == z and not _is_nitro_or_nitroso(a)}
    runs, seen = [], set()
    for start in sorted(atoms):
        if start in seen:
            continue
        component, stack = set(), [start]
        while stack:
            a = stack.pop()
            if a in component:
                continue
            component.add(a)
            stack.extend(n for n in graph[a] if n in atoms)
        seen |= component
        ends = [a for a in component if sum(n in component for n in graph[a]) <= 1]
        if len(component) > 1 and (len(ends) != 2 or any(sum(n in component for n in graph[a]) > 2 for a in component)):
            return None
        run, previous = [ends[0]], None
        while True:
            following = [n for n in graph[run[-1]] if n in component and n != previous]
            if not following:
                break
            previous = run[-1]
            run.append(following[0])
        if any(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0 for a, b in zip(run, run[1:])):
            return None
        runs.append(run)
    if len(runs) < 2:
        return None
    best = max(runs, key=lambda run: (len(run), sum(len(graph[a]) - sum(n in run for n in graph[a]) for a in run)))
    if not _joined_by_chalcogens(mol, graph, best, runs):
        return None
    return z, best


def _joined_by_chalcogens(mol, graph, best, runs):
    """Whether every other run is reached from `best` only through O, S, Se or Te atoms (which rank below the runs)."""
    others = {a for run in runs if run is not best for a in run}
    parent = {a: None for a in best}
    queue = list(best)
    while queue:
        a = queue.pop(0)
        for n in graph[a]:
            if n not in parent:
                parent[n] = a
                queue.append(n)
    for run in runs:
        if run is best:
            continue
        a = parent.get(run[0])
        while a is not None and a not in best:
            if a not in others and mol.GetAtomWithIdx(a).GetAtomicNum() not in (8, 16, 34, 52):
                return False
            a = parent[a]
    return True


def _unsaturated_parent(z, size, ene, yne, grouped):
    """'disilyne', 'triazene', 'triaz-1-ene', 'pentaaza-1,3-diene' (P-14.3.4.2(d), P-31.1.4): the 'ane' of the chain hydride
    becomes 'ene' or 'yne'; the locants are omitted for a dinuclear chain and for an unsubstituted monounsaturated one."""
    root = multiplying_prefix(size) + _STEMS[z][:-3]
    unsubstituted = not grouped
    if size == 2 or (size == 3 and len(ene) + len(yne) == 1 and unsubstituted):
        return root + ("yne" if yne and not ene else "ene" if ene and not yne else _unsat_words(ene, yne)), True
    body, needs_a = unsaturation_suffix(ene, yne)
    return f"{root}{'a' if needs_a else ''}-{body}", False


def _unsat_words(ene, yne):
    from ._common import unsaturation_suffix

    return unsaturation_suffix(ene, yne)[0]


def name_hydride_chain(mol, graph, halogens, aromatic_atoms):
    found = _chain_atoms(mol, graph)
    unsaturated = False
    if found is None:
        found = _chain_atoms(mol, graph, allow_double=True, allow_triple=True)
        unsaturated = found is not None
    alternating = None
    if found is None:
        alternating = _alternating_chain(mol, graph)
        if alternating is None:
            found = _split_chain(mol, graph)
            if found is None:
                return None
            z, chain = found
        else:
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
        ene = [i + 1 for i in range(len(candidate) - 1) if mol.GetBondBetweenAtoms(candidate[i], candidate[i + 1]).GetBondTypeAsDouble() == 2.0]
        yne = [i + 1 for i in range(len(candidate) - 1) if mol.GetBondBetweenAtoms(candidate[i], candidate[i + 1]).GetBondTypeAsDouble() == 3.0]
        lam = {i + 1: n for i, atom in enumerate(candidate) if (n := nonstandard_bonding(mol.GetAtomWithIdx(atom)))}
        omit = not lam and (
            len(chain) == 1
            or (len(chain) == 2 and sum(len(info["locants"]) for info in grouped.values()) == 1)
            or all(mol.GetAtomWithIdx(a).GetTotalNumHs() == 0 for a in chain)
        )
        if alternating is not None:
            parent = _alternating_parent(alternating[0], alternating[1], (len(chain) + 1) // 2, ene, yne)
            cited = sum(len(info["locants"]) for info in grouped.values())
            fully_substituted = len(grouped) == 1 and all(mol.GetAtomWithIdx(a).GetTotalNumHs() == 0 for a in chain)
            omit = not lam and ((len(chain) == 3 and cited == 1) or fully_substituted)
        elif unsaturated:
            parent, omit = _unsaturated_parent(z, len(chain), ene, yne, grouped)
        else:
            parent = _STEMS[z] if len(chain) == 1 else f"{multiplying_prefix(len(chain))}{_STEMS[z]}"
            if z == 7 and len(chain) == 2:
                parent = "hydrazine"
        lam_text = ",".join(f"{p}\u03bb{lam[p]}" for p in sorted(lam))
        if lam_text:
            parent = f"{lam_text}-{parent}"
        if omit and len(grouped) > 1:
            entries = [(n, info["compound"]) for n, info in grouped.items() for _ in info["locants"]]
            name = format_mononuclear_prefixes(entries) + parent
        else:
            prefixes = format_substituent_prefixes(grouped, omit_locants=omit)
            name = prefixes + ("-" if prefixes and lam_text else "") + parent
        key = (sorted(lam), [-lam[p] for p in sorted(lam)], ene, yne, locant_set, citation, name)
        if best is None or key < best[0]:
            best = (key, name)
    return best[1]
