"""Anionic centers on an unbranched chain of one repeated heteroatom (P-72.2.2.1,
P-72.2.2.2.4, P-72.6.3): hydrazine-1,1-diide, trioxidanide, disulfanidyl, ... The 'ide'
centers take the lowest locants, then the free valence and the substituent prefixes
(P-31.1.4); locants that cannot be ambiguous are omitted (P-14.3.4).
"""

from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents, multiplied_word
from ._substituents import format_substituent_prefixes, name_branch

_CHAIN_STEMS = {
    7: ("azane", "hydrazine"),
    8: ("oxidane", None),
    14: ("silane", None),
    15: ("phosphane", None),
    16: ("sulfane", None),
    32: ("germane", None),
    33: ("arsane", None),
    34: ("selane", None),
    50: ("stannane", None),
    52: ("tellane", None),
}
_CHALCOGENS = {8, 16, 34, 52}
_STANDARD = {7: 3, 8: 2, 14: 4, 15: 3, 16: 2, 32: 4, 33: 3, 34: 2, 50: 4, 52: 2}
_ENDING = {1.0: "yl", 2.0: "ylidene", 3.0: "ylidyne"}


def is_chain_anion(mol):
    for atom in mol.GetAtoms():
        if atom.HasProp("_anion_word") and any(n.GetAtomicNum() == atom.GetAtomicNum() for n in atom.GetNeighbors()):
            return atom.GetAtomicNum() in _CHAIN_STEMS
    return False


def name_chain_anion(neutral):
    centers = [a for a in neutral.GetAtoms() if a.HasProp("_anion_word")]
    z = centers[0].GetAtomicNum()
    if any(a.GetAtomicNum() != z for a in centers) or z not in _CHAIN_STEMS:
        raise UnsupportedStructure("the anionic centers of a heteroatom chain must be one element")
    chain = _chain_atoms(neutral, centers[0].GetIdx(), z, None)
    if any(a.GetIdx() not in chain for a in centers):
        raise UnsupportedStructure("every anionic center must lie on the one heteroatom chain")
    return _compose(neutral, chain, centers, z, None)


def chain_prefix(neutral, root, coming_from):
    """'disulfanidyl'-type prefix for a heteroatom chain attached through `root`, or None."""
    z = neutral.GetAtomWithIdx(root).GetAtomicNum()
    if z not in _CHAIN_STEMS:
        return None
    chain = _chain_atoms(neutral, root, z, coming_from)
    centers = [neutral.GetAtomWithIdx(a) for a in chain if neutral.GetAtomWithIdx(a).HasProp("_anion_word")]
    if len(chain) < 2 or not centers:
        return None
    if root not in (chain[0], chain[-1]):
        raise UnsupportedStructure("an anionic heteroatom chain must be attached through an end atom")
    text = _compose(neutral, chain, centers, z, (root, coming_from))
    return text, any(ch.isdigit() for ch in text) or "(" in text


def _compose(neutral, chain, centers, z, free):
    words = {a.GetProp("_anion_word") for a in centers}
    if len(words) != 1:
        raise UnsupportedStructure("ide and uide centers in one chain are not supported yet")
    word = words.pop()
    graph = adjacency(neutral)
    halogens = halogen_substituents(neutral)
    aromatic_atoms = frozenset(a.GetIdx() for a in neutral.GetAtoms() if a.GetIsAromatic())
    chain_set = set(chain)
    blocked = free[1] if free else None
    best = None
    for ordered in (chain, chain[::-1]):
        position = {atom: i + 1 for i, atom in enumerate(ordered)}
        entries = {}
        for atom in ordered:
            for n in graph[atom]:
                if n in chain_set or n == blocked:
                    continue
                entries.setdefault(position[atom], []).append(
                    name_branch(graph, n, atom, halogens, aromatic_atoms, mol=neutral, unsaturated=True)
                )
        locants = sorted(position[a.GetIdx()] for a in centers for _ in range(int(a.GetProp("_anion_charge"))))
        free_locant = position[free[0]] if free else 0
        grouped = group_substituents(entries)
        key = (
            tuple(locants),
            free_locant,
            sorted(loc for info in grouped.values() for loc in info["locants"]),
            sorted(grouped),
        )
        if best is None or key < best[0]:
            best = (key, locants, grouped, position, free_locant)
    _, locants, grouped, position, free_locant = best
    lambdas = {position[a.GetIdx()]: a.GetProp("_anion_lambda") for a in centers if a.HasProp("_anion_lambda")}
    stem, retained = _CHAIN_STEMS[z]
    omit = z in _CHALCOGENS and not lambdas
    prefixes = format_substituent_prefixes(
        grouped, omit_locants=omit or (not lambdas and _one_arrangement_of(grouped, locants, free_locant, len(chain), z, neutral, chain, position))
    )
    n = len(chain)
    parent = retained if (retained and n == 2) else multiplied_word(n, stem)
    if lambdas:
        parent = ",".join(f"{loc}λ{lambdas[loc]}" for loc in sorted(lambdas)) + "-" + parent
    count_word = multiplied_word(len(locants), word)
    if free:
        count_word = count_word[:-1]
    if count_word[0] in "aeiouy":
        parent = parent[:-1]
    locant_text = "" if omit else f"-{','.join(str(x) for x in locants)}-"
    text = f"{prefixes}{parent}{locant_text}{count_word}"
    if free:
        order = neutral.GetBondBetweenAtoms(*free).GetBondTypeAsDouble()
        text += ("" if omit else f"-{free_locant}-") + _ENDING[order]
    return text


def _chain_atoms(mol, start, z, blocked):
    seen = [start]
    frontier = [start]
    while frontier:
        current = frontier.pop()
        for n in mol.GetAtomWithIdx(current).GetNeighbors():
            if n.GetAtomicNum() == z and n.GetIdx() not in seen and n.GetIdx() != blocked:
                seen.append(n.GetIdx())
                frontier.append(n.GetIdx())
    degrees = {a: sum(1 for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in seen) for a in seen}
    ends = [a for a, d in degrees.items() if d <= 1]
    if any(d > 2 for d in degrees.values()) or (len(seen) > 1 and len(ends) != 2) or any(
        mol.GetAtomWithIdx(a).IsInRing() for a in seen
    ):
        raise UnsupportedStructure("a branched or cyclic heteroatom chain is not supported here")
    if len(seen) == 1:
        return seen
    ordered = [ends[0]]
    while len(ordered) < len(seen):
        nxt = next(
            n.GetIdx()
            for n in mol.GetAtomWithIdx(ordered[-1]).GetNeighbors()
            if n.GetIdx() in seen and n.GetIdx() not in ordered
        )
        ordered.append(nxt)
    return ordered


def _one_arrangement_of(grouped, locants, free_locant, size, z, mol, chain, position):
    """P-14.3.4.4: whether moving a prefix to another chain atom, or exchanging two of them, can only give the same
    compound once the cited anionic centers and free valence hold their atoms, so that the locants of the prefixes
    are omitted ('acetylhydrazine-1,1-diide')."""
    from ._hydride_chain import _one_arrangement

    if any(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0 for a, b in zip(chain, chain[1:])):
        return False
    entries = {}
    for name, info in grouped.items():
        for locant in info["locants"]:
            entries.setdefault(locant, []).append((name, info["compound"]))
    taken = [*locants, *([free_locant] if free_locant else [])]
    capacity = [
        _STANDARD[z] - sum(1 for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in chain) - taken.count(position[a])
        for a in sorted(chain, key=position.get)
    ]
    return len(capacity) == size and _one_arrangement(entries, capacity)
