"""Naming of compound (branched) substituent groups, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-29.4.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): "A
  compound substituted substituent group is formed by substituting one or
  more simple substituents into another simple substituent that is
  considered as the principal chain." The free valence (the atom attached to
  the parent chain or ring) is always locant 1 of that principal chain;
  remaining branches are cited as nested substituent prefixes, with identical
  branches grouped under a multiplying prefix.
- P-46 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): the
  principal chain of a substituent group is chosen by the same criteria as a
  parent hydride's (P-44.3/P-45.2) — longest chain, then most substituents,
  then lowest locant set, then lowest locants in citation order — except that
  locant 1 is fixed at the free valence, so there is no choice of numbering
  direction the way there is for a parent hydride.
- P-14.2.2 (Chapter P-1): 'bis', 'tris', 'tetrakis', ... multiply identical
  compound substituent prefixes, instead of 'di', 'tri', 'tetra', ..., to
  avoid ambiguity with a substituent's own internal multiplying prefixes.
- P-14.5.2 (Chapter P-1): alphanumerical order is based on a substituent
  prefix's complete name, so a compound substituent like '(1-methylpropyl)'
  alphabetizes under 'm' (from 'methylpropyl'), ignoring the enclosing
  parentheses and locants. Halogeno prefixes (P-35.2.1, see `_common.py`)
  alphabetize the same way, under their own name ('bromo', 'chloro',
  'fluoro', 'iodo') — no special-casing needed.
- P-35.2.1 (Chapter P-3): a halogen atom (F, Cl, Br, I) directly attached to
  a chain/ring atom is itself a simple substituent group ('fluoro', 'chloro',
  'bromo', 'iodo') with no locants or nested prefixes of its own — the base
  case in `name_branch` below.

Cyclic substituent groups (P-29.3.3) are out of scope and raise
UnsupportedStructure, except the minimal case added for
`tasks/ring-substituent-chain-suffix.md` (2026-08-25): a plain, unsubstituted
saturated monocyclic ring hanging off the parent chain (e.g. "cyclohexyl" in
cyclohexylmethanol) is recognized by `_simple_ring_substituent` and named
directly ("cyclo" + `alkyl_name`), without walking into
`_longest_chains_from_root`'s cycle-detection rejection. A ring bearing its
own substituent, an unsaturated ring, or a polycyclic/spiro ring as a
substituent all remain out of scope and still raise `UnsupportedStructure`
via that same cycle-detection path (see `_simple_ring_substituent`'s own
docstring for exactly which shapes it recognizes).
"""

import re

from ._common import UnsupportedStructure, lowest_locant_set
from ._numerals import alkyl_name, multiplying_prefix

_LEADING_LOCANTS_RE = re.compile(r"^[\d,\-]+")


def alpha_sort_key(name: str) -> str:
    return _LEADING_LOCANTS_RE.sub("", name).lower()


def format_substituent_prefixes(grouped, omit_locants: bool = False) -> str:
    """grouped: {name -> {"locants": [int, ...], "compound": bool}}. Return
    the assembled, alphanumerically ordered prefix string (P-14.5.2), ready to
    prepend to a parent name; '' if grouped is empty.

    `omit_locants`: for a mononuclear parent hydride (P-14.3.4.2(a)), every
    substituent's locant is always '1' and never cited, no matter how many
    substituents there are — unlike the ordinary case, where a locant is
    droppable only when every substituent's is unambiguous without it."""
    parts = []
    for name in sorted(grouped, key=alpha_sort_key):
        info = grouped[name]
        locants = sorted(info["locants"])
        multiplier = multiplying_prefix(len(locants), compound=info["compound"]) if len(locants) > 1 else ""
        display_name = f"({name})" if info["compound"] else name
        if omit_locants:
            parts.append(f"{multiplier}{display_name}")
        else:
            loc_str = ",".join(str(loc) for loc in locants)
            parts.append(f"{loc_str}-{multiplier}{display_name}")
    return "".join(parts) if omit_locants else "-".join(parts)


def _group_substituents(entries):
    """entries: iterable of (locant, name, is_compound)."""
    grouped = {}
    for locant, name, is_compound in entries:
        info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
        info["locants"].append(locant)
    return grouped


def _longest_chains_from_root(graph, root, coming_from, halogens):
    """All maximum-length simple paths starting at `root`, extending into the
    subtree away from `coming_from` (P-46: the free valence is fixed at
    locant 1, so only one direction of travel is possible, unlike a parent
    hydride's chain). `halogens` atoms are excluded from the walk (P-35.2.1:
    a halogen is a terminal substituent, never a chain-extending atom), the
    same way `carbon_adjacency` excludes them from a parent hydride's chain
    search."""
    best_length = 0
    best_paths = []

    def walk(node, previous, path, path_set):
        nonlocal best_length, best_paths
        neighbors = [n for n in graph[node] if n != previous and n not in halogens]
        if not neighbors:
            if len(path) > best_length:
                best_length, best_paths = len(path), [list(path)]
            elif len(path) == best_length:
                best_paths.append(list(path))
            return
        for neighbor in neighbors:
            if neighbor in path_set:
                raise UnsupportedStructure(
                    "cyclic substituent groups are not supported yet (see "
                    "P-29.3.3, P-46 for cyclic substituent groups)"
                )
            path.append(neighbor)
            path_set.add(neighbor)
            walk(neighbor, node, path, path_set)
            path_set.discard(neighbor)
            path.pop()

    walk(root, coming_from, [root], {root})
    return best_paths


def _candidate_key(grouped):
    """Sort key implementing P-46's analogue of P-45.2.1-P-45.2.3, most
    preferred first (chain length is fixed by the caller, so it is not part
    of this key)."""
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    return -total_count, locant_set, citation_locants


def _simple_ring_substituent(graph, root, coming_from):
    """If the branch hanging off `root` (away from `coming_from`) is a
    single, simple, unsubstituted saturated monocyclic ring with `root` as
    its only attachment point, return the ring size; else None (a
    non-ring branch, a ring bearing its own substituent, or any
    polycyclic/spiro/fused shape all fall through to the ordinary
    chain-walk in `name_branch`, which raises `UnsupportedStructure` via
    `_longest_chains_from_root`'s cycle-detection check)."""
    ring_neighbors = [n for n in graph[root] if n != coming_from]
    if len(ring_neighbors) != 2:
        return None
    visited = {root}
    previous, current = root, ring_neighbors[0]
    while current != root:
        if current in visited:
            return None
        visited.add(current)
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) != 1:
            return None
        previous, current = current, neighbors[0]
    return len(visited)


def name_branch(graph, root, coming_from, halogens=None):
    """Name the substituent group hanging off `root`, reached from
    `coming_from` (the parent chain/ring atom). Returns (name, is_compound);
    is_compound is True iff the name carries its own locants/nested prefixes
    (P-29.4) and should be parenthesized when cited as a prefix.

    `halogens`: {atom_idx -> prefix name} (see `_common.halogen_substituents`).
    If `root` is itself a halogen atom, it's a simple substituent with no
    locants or nested prefixes of its own (P-35.2.1) — returned directly,
    with no recursion."""
    halogens = halogens or {}
    if root in halogens:
        return halogens[root], False

    ring_size = _simple_ring_substituent(graph, root, coming_from)
    if ring_size is not None:
        return "cyclo" + alkyl_name(ring_size), False

    chains = _longest_chains_from_root(graph, root, coming_from, halogens)
    chain_length = len(chains[0])

    best_key = None
    best_name = None
    best_compound = None
    for chain in chains:
        chain_set = set(chain)
        entries = []
        for position, atom in enumerate(chain, start=1):
            previous = chain[position - 2] if position > 1 else coming_from
            for branch_root in graph[atom]:
                if branch_root == previous or branch_root in chain_set:
                    continue
                sub_name, sub_compound = name_branch(graph, branch_root, atom, halogens)
                entries.append((position, sub_name, sub_compound))
        grouped = _group_substituents(entries)
        if grouped and chain_length == 1:
            # P-14.3.4.2(a): the branch's own chain is a single (mononuclear)
            # atom, so any substituent on it has no other possible position
            # and its locant is never cited, e.g. '(chloromethyl)', not
            # '(1-chloromethyl)'.
            name, is_compound = format_substituent_prefixes(grouped, omit_locants=True) + alkyl_name(1), True
        elif grouped:
            name, is_compound = format_substituent_prefixes(grouped) + alkyl_name(chain_length), True
        else:
            name, is_compound = alkyl_name(chain_length), False
        key = _candidate_key(grouped) + (name,)
        if best_key is None or key < best_key:
            best_key, best_name, best_compound = key, name, is_compound

    return best_name, best_compound
