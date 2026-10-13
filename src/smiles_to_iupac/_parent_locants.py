"""Numberings of an aldehyde or ketone parent that its name expresses (P-31.1.4, P-44.1.2.2, P-44.3): the
atom-to-locant maps of the principal carbonyl carbons, including every tied numbering."""

import re

from ._common import UnsupportedStructure, adjacency
from ._diester_ring_diyl import _system_of, evaluate_skeleton

_MAX_AUTOMORPHISMS = 50000
_LOCANT = re.compile(r"(\d+)([a-z]*)([\'′]*)")
_CHAIN_GROUP = re.compile(
    r"(?P<stem>.+)-(?P<locants>\d+(?:,\d+)*)-(?P<multiplier>(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca|undeca|dodeca)?)yl"
)
_PHENYLENE = re.compile(r"(?P<prefixes>.*-)?(?P<locants>\d+(?:,\d+)*)-phenylene")


def locant_key(locant):
    match = _LOCANT.fullmatch(str(locant))
    if match is None:
        raise UnsupportedStructure(f"the locant {locant!r} cannot be ordered")
    return int(match.group(1)), match.group(2), len(match.group(3))


def parent_numberings(parent, principal, expected_name):
    """[{carbon: locant}] for the principal carbonyl carbons of `parent` under every numbering that gives the
    systematic name `expected_name`; empty when no ring system or chain reproduces it."""
    graph = adjacency(parent)
    aldehyde = parent.GetAtomWithIdx(next(iter(principal))).GetTotalNumHs() > 0
    found = _ring_numberings(parent, graph, principal, aldehyde, expected_name)
    found += _chain_numberings(parent, graph, principal, aldehyde, expected_name)
    return _with_symmetry(parent, found)


def _oxygen(parent, carbon):
    return next(n.GetIdx() for n in parent.GetAtomWithIdx(carbon).GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1)


def _ring_numberings(parent, graph, principal, aldehyde, expected):
    anchors = {}
    for carbon in principal:
        atom = parent.GetAtomWithIdx(carbon)
        if atom.IsInRing():
            anchors[carbon] = carbon
        elif aldehyde:
            ring_neighbours = [n.GetIdx() for n in atom.GetNeighbors() if n.IsInRing()]
            if len(ring_neighbours) == 1:
                anchors[carbon] = ring_neighbours[0]
    systems = {}
    for carbon, anchor in anchors.items():
        rings, atoms = _system_of(parent, anchor)
        systems.setdefault(frozenset(atoms), (rings, atoms, []))[2].append(carbon)
    found = []
    for rings, atoms, carbons in systems.values():
        attach = [anchors[c] for c in carbons]
        blocked = {_oxygen(parent, c) for c in carbons} | {c for c in carbons if anchors[c] != c}
        suffix = "carbaldehyde" if aldehyde else "one"
        try:
            named = evaluate_skeleton(parent, graph, "ring", rings, atoms, attach, blocked, suffix)
        except UnsupportedStructure:
            continue
        if named is not None and _ring_name_is(named[1], expected, suffix, len(carbons)):
            found.append({c: named[2][anchors[c]] for c in carbons})
    return found


def _ring_name_is(group_name, expected, suffix, count):
    if group_name == expected:
        return True
    match = _PHENYLENE.fullmatch(group_name)
    if match is None or suffix != "carbaldehyde":
        return False
    from ._numerals import multiplying_prefix

    prefixes = (match.group("prefixes") or "")[:-1]
    multiplier = multiplying_prefix(count) if count > 1 else ""
    return expected == f"{prefixes}benzene-{match.group('locants')}-{multiplier}carbaldehyde"


def _chain_numberings(parent, graph, principal, aldehyde, expected):
    acyclic = {a.GetIdx() for a in parent.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
    attach = [c for c in principal if c in acyclic]
    if not attach:
        return []
    forest = {a: [n for n in graph[a] if n in acyclic] for a in acyclic}
    blocked = {_oxygen(parent, c) for c in attach}
    found = []
    for path in _end_to_end_paths(forest, attach):
        try:
            named = evaluate_skeleton(parent, graph, "chain", [path], set(path), attach, blocked, "one")
        except UnsupportedStructure:
            continue
        if named is not None and _chain_name_is(named[1], expected, aldehyde):
            found.append({c: named[2][c] for c in attach})
    return found


def _end_to_end_paths(forest, attach):
    """Carbon chains between two chain ends (carbons with at most one carbon neighbour) that hold every `attach`."""
    ends = [a for a, neighbours in forest.items() if len(neighbours) <= 1]
    required = set(attach)
    paths = []

    def extend(path):
        last = path[-1]
        closed = len(forest[last]) <= 1 and (path[0] < last or len(path) == 1 and not forest[last])
        if closed and required <= set(path):
            paths.append(list(path))
        for n in forest[last]:
            if n not in path:
                extend(path + [n])

    for end in ends:
        extend([end])
    return paths


def _chain_name_is(group_name, expected, aldehyde):
    match = _CHAIN_GROUP.fullmatch(group_name)
    if match is None:
        return aldehyde and group_name.endswith("yl") and expected == group_name[:-2] + "anal"
    stem, multiplier = match.group("stem"), match.group("multiplier")
    if aldehyde:
        return expected == f"{stem}{multiplier}al"
    return expected == f"{stem}-{match.group('locants')}-{multiplier}one"


def _with_symmetry(parent, numberings):
    """Every numbering that the automorphisms of the parent carry the found ones into: tied numberings give the same
    parent name, so the derivative groups decide between them."""
    matches = parent.GetSubstructMatches(parent, uniquify=False, useChirality=True, maxMatches=_MAX_AUTOMORPHISMS)
    seen = set()
    result = []
    for base in numberings:
        for match in matches:
            image = {match[atom]: locant for atom, locant in base.items()}
            key = frozenset(image.items())
            if key not in seen:
                seen.add(key)
                result.append(image)
    return result
