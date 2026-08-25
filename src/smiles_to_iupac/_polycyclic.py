"""Generic von Baeyer naming engine for saturated polycyclic hydrocarbons
whose skeleton reduces to exactly `2*(ring_count-1)` skeletal atoms of degree
3 (and none of higher degree), joined by `3*(ring_count-1)` bridges: a main
bicyclic system (P-23.2.2/P-23.2.3, two main bridgeheads and three bridges
between them) plus `ring_count-2` further *independent* secondary bridges
among the remaining branch atoms, each attaching only to atoms already placed
by the main system (P-23.2.5.1's "independent secondary bridge" notion,
generalized to more than one such bridge exactly as VB-6/VB-7 do -- see
below).

This is the consolidated form of what were previously separate, near-
identical `_tricyclic.py` (ring_count=3, one secondary bridge) and
`_tetracyclic.py` (ring_count=4, two secondary bridges) modules: both
`find_*_core` functions and both `name_*cycloalkane` functions differed only
in a hardcoded secondary-bridge count, so this module takes that count as the
`ring_count` parameter instead. `_bicyclic.py` (ring_count=2, zero secondary
bridges, so no branch-atom "hopping" is ever needed) is intentionally *not*
folded in here: its two-bridgehead representation and simpler search loop
are different enough in shape that unifying them would obscure both, per
`tasks/von-baeyer-engine-consolidation.md`'s decision to consolidate only
once a third near-identical module (this one's pentacyclic addition) made
the duplication cost concrete.

Per the IUPAC 2013 Recommendations ("the Blue Book", Chapter P-2,
https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) and W.H. Powell, "Extension
and Revision of the von Baeyer System for Naming Polycyclic Compounds", Pure
Appl. Chem. 71 (1999) 513-529 (mirrored at
https://iupac.qmul.ac.uk/vonBaeyer/vb6.html and .../vb7.html):

- P-23.2.1 (max main ring) and P-23.2.4 (max main bridge): the main ring and
  main bridge are chosen to include as many skeletal atoms as possible,
  taking priority over every tie-break below.
- P-23.2.6.2.1 (symmetry): remaining ties in the main-ring/main-bridge choice
  are broken by dividing the main ring as symmetrically as possible between
  its two segments.
- VB-6 ("Naming Hydrocarbon Polycyclic Systems"): the independent secondary
  bridges are chosen, among the bridges left over once the main system is
  fixed, to maximize the largest one's length, then the next, and so on; they
  are *cited* in that same decreasing-size order (equal-length ties cited in
  ascending locant order -- see the tetracyclo[4.4.2.2^2,5.2^7,10]hexadecane
  worked example, where the two equal-length secondary bridges are cited
  (2,5) before (7,10)).
- VB-7 ("Numbering of Secondary Bridges"): after the main ring and main
  bridge are numbered, the independent secondary bridges are numbered in
  order of the highest-numbered bridgehead each attaches to (higher first),
  ties broken by longer bridge first, then lower attachment locant --
  independent of, and possibly different from, citation order above. Within
  a single secondary bridge, numbering starts adjacent to its *higher*-
  numbered attachment point, but the pair is still cited lower-locant-first.
- P-23.2.6.2.4/.2.5-style tie-breaks: once main system and secondary-bridge
  lengths are fixed, prefer the decomposition with the lowest combined
  locant set for all secondary-bridge attachment points, then the lowest
  locants in citation order.
- P-14.4 / P-45.2 (Chapter P-1 / P-4): any choice still left open is settled
  by lowest locants to substituents, as in every other ring module here.

Only *independent* secondary bridges are supported (each one's own two
attachment points must both already lie on the main ring/main bridge) --
"dependent" secondary bridges, whose attachment point lies on another
secondary bridge, are not detected by the search below and simply yield no
valid decomposition. Propellane-like topologies (any branch atom of degree 4
or higher) are also out of scope; see `_tricyclic.py`'s
`find_propellane_core`/`name_propellane` for the one degree-4 case that is
supported (tricyclic only). Any topology that doesn't reduce to the shape
above raises `UnsupportedStructure` by simply finding no valid decomposition.

Flagship validation cases:
- tricyclic (ring_count=3): adamantane, tricyclo[3.3.1.1^3,7]decane; see
  `_tricyclic.py`.
- tetracyclic (ring_count=4): quadricyclane, tetracyclo[3.2.0.0^2,7.0^4,6]-
  heptane (C7H8; PubChem CID 78961 gives canonical SMILES C1C2C3C2C4C1C34 and
  this exact PubChem-computed IUPAC name -- cross-checked externally via the
  PubChem PUG REST API, not from memory, and independently re-verified via
  RDKit's reported formula/degree sequence).
- pentacyclic (ring_count=5): cubane, pentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]-
  octane (C8H8, eight branch atoms all of degree 3, cyclomatic number 5;
  cross-checked against PubChem CID 136090's computed IUPAC name via the PUG
  REST API).
"""

from itertools import combinations, permutations

from ._cyclic import _group, _substituents_for_ring
from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes


def _strip_leaves(graph):
    """Repeatedly remove degree-<=1 atoms; what remains is the 2-connected
    ring-system core (see `_bicyclic.py`'s identical helper)."""
    core = {atom: set(neighbors) for atom, neighbors in graph.items()}
    changed = True
    while changed:
        changed = False
        for leaf in [atom for atom, neighbors in core.items() if len(neighbors) <= 1]:
            for neighbor in core[leaf]:
                core[neighbor].discard(leaf)
            del core[leaf]
            changed = True
    return core


def _walk_to_branch(core, branch_atoms, start, first):
    """Walk from `start` via its neighbor `first` until another branch atom
    (one of `branch_atoms`) is reached, returning (branch_atom, internal_path)
    with the path ordered nearest-`start`-first, or None if the walk loops
    back to `start` itself without reaching a different branch atom (an
    invalid/self-looping bridge)."""
    if first in branch_atoms:
        return first, []
    path = [first]
    previous, current = start, first
    while True:
        neighbors = core[current] - {previous}
        if len(neighbors) != 1:
            return None
        next_atom = next(iter(neighbors))
        if next_atom == start:
            return None
        if next_atom in branch_atoms:
            return next_atom, path
        path.append(next_atom)
        previous, current = current, next_atom
    return None


def find_polycyclic_core(mol, ring_count):
    """Return (branch_atoms, bridges) if `mol`'s carbon skeleton, after
    stripping acyclic branches, reduces to exactly `2*(ring_count-1)` branch
    atoms of degree 3 (and no atom of higher degree) joined by
    `3*(ring_count-1)` bridges (the graph's cyclomatic number is
    `ring_count`) -- else None.

    `bridges` is a list of `(u, v, path)` tuples, one per bridge, where
    `path` is the list of internal atom indices ordered nearest-`u`-first.
    Branch-atom pairs need not each have exactly one bridge between them:
    doubled bridges between the same pair are a real, connected topology and
    are detected here the same way.

    Mirrors `_bicyclic.find_bicyclic_core`'s leaf-stripping + cyclomatic-
    number approach rather than trusting RDKit's `NumRings()` (see that
    function's docstring for why)."""
    graph = adjacency(mol)
    core = _strip_leaves(graph)
    if not core:
        return None
    vertices = len(core)
    edge_count = sum(len(neighbors) for neighbors in core.values()) // 2
    if edge_count - vertices + 1 != ring_count:
        return None
    if any(len(neighbors) not in (2, 3) for neighbors in core.values()):
        return None
    branch_atoms = {atom for atom, neighbors in core.items() if len(neighbors) == 3}
    if len(branch_atoms) != 2 * (ring_count - 1):
        return None

    half_walks = []
    for u in branch_atoms:
        for first in core[u]:
            result = _walk_to_branch(core, branch_atoms, u, first)
            if result is None:
                return None
            v, path = result
            half_walks.append((u, v, path))
    if len(half_walks) != 6 * (ring_count - 1):
        return None

    bridges = []
    used = [False] * len(half_walks)
    for i, (u, v, path) in enumerate(half_walks):
        if used[i]:
            continue
        reversed_path = list(reversed(path))
        match = None
        for j in range(i + 1, len(half_walks)):
            if used[j]:
                continue
            u2, v2, path2 = half_walks[j]
            if u2 == v and v2 == u and path2 == reversed_path:
                match = j
                break
        if match is None:
            return None
        used[i] = used[match] = True
        bridges.append((u, v, path))
    if len(bridges) != 3 * (ring_count - 1):
        return None

    return branch_atoms, bridges


def _adjacency_by_atom(branch_atoms, bridges):
    """{atom: [(bridge_index, other_atom, path_from_atom_to_other), ...]}."""
    adj = {atom: [] for atom in branch_atoms}
    for idx, (u, v, path) in enumerate(bridges):
        adj[u].append((idx, v, path))
        adj[v].append((idx, u, list(reversed(path))))
    return adj


def _composite_bridges(adj, start, end, others):
    """Yield (bridge_indices_used, atom_path) for every simple walk from
    `start` to `end` that uses any subset (in any order, zero or more) of
    `others` as intermediate branch atoms, each hop a real bridge with no
    bridge index reused and no branch atom visited twice. `atom_path` is the
    full ordered list of atoms strictly between `start` and `end`, including
    any intermediate branch atom(s) -- so when a composite is later used as a
    main-system bridge, its intermediate branch atom(s) fall naturally into
    the resulting numbering.

    A bounded-depth backtracking walk over an arbitrary number of "other"
    branch atoms, so it works unchanged for any `ring_count`."""
    others = set(others)

    def walk(current, used_idxs, used_atoms, path):
        for idx, v, seg in adj[current]:
            if idx in used_idxs:
                continue
            if v == end:
                yield used_idxs + (idx,), path + seg
            elif v in others and v not in used_atoms:
                yield from walk(v, used_idxs + (idx,), used_atoms | {v}, path + seg + [v])

    yield from walk(start, (), frozenset(), [])


def iter_polycyclic_candidates(core, ring_count):
    """Yield `(full_order, parent, outer_key)` for every von Baeyer-valid
    numbering of `core` (P-23.2.1-P-23.2.6, VB-6/VB-7) -- every choice of
    main bridgeheads/main-system decomposition and, among candidates tied on
    ring/bridge shape, every remaining numbering degree of freedom (main-ring
    direction, order among equal-length bridges).

    `outer_key` ranks the topology/descriptor choice (bridge-length shape,
    secondary-bridge lengths and their locants) -- everything that must be
    decided *before* substituent (or, as of
    `tasks/heteroatom-skeleton-expansion.md`, 2026-08-26, a skeletal
    heteroatom's) locants can break a remaining tie, since `parent` itself
    (unlike `_bicyclic.py`'s fixed 'bicyclo[x.y.z]alkane') already encodes a
    choice that must be settled first. A caller combines
    `outer_key + _candidate_key(parent, substituents, ...)` to get the full
    ranking key, mirroring `_bicyclic.iter_bicyclic_numberings`'s split
    (shared with `_von_baeyer_heteroatom.py` for the same reason)."""
    branch_atoms, bridges = core
    adj = _adjacency_by_atom(branch_atoms, bridges)
    secondary_count = ring_count - 2
    ring_prefix = numerical_term(ring_count) + "cyclo"

    for bh1, bh2 in combinations(branch_atoms, 2):
        others = tuple(a for a in branch_atoms if a not in (bh1, bh2))
        for start, other in ((bh1, bh2), (bh2, bh1)):
            composites = list(_composite_bridges(adj, start, other, others))
            for i, j, k in combinations(range(len(composites)), 3):
                idxs = [composites[i][0], composites[j][0], composites[k][0]]
                paths = [composites[i][1], composites[j][1], composites[k][1]]
                used_idxs = set(idxs[0]) | set(idxs[1]) | set(idxs[2])
                if len(used_idxs) != len(idxs[0]) + len(idxs[1]) + len(idxs[2]):
                    continue
                remaining = set(range(len(bridges))) - used_idxs
                if len(remaining) != secondary_count:
                    continue

                a, b, c = sorted((len(p) for p in paths), reverse=True)
                # P-23.2.1 (max main ring) > P-23.2.4 (max main bridge) >
                # P-23.2.6.2.1 (main ring divided as symmetrically as possible).
                shape_key = (-(a + b), -c, a - b)

                for perm in permutations(paths):
                    if not (len(perm[0]) >= len(perm[1]) >= len(perm[2])):
                        continue
                    main_ring_first, main_ring_second, main_bridge = perm
                    main_order = (
                        [start] + main_ring_first + [other]
                        + list(reversed(main_ring_second)) + main_bridge
                    )
                    position = {atom: idx + 1 for idx, atom in enumerate(main_order)}

                    secondary = []
                    all_independent = True
                    for sec_idx in remaining:
                        sec_u, sec_v, sec_path_uv = bridges[sec_idx]
                        if sec_u not in position or sec_v not in position:
                            all_independent = False
                            break
                        x, y = position[sec_u], position[sec_v]
                        lo, hi = min(x, y), max(x, y)
                        hi_atom = sec_u if x > y else sec_v
                        secondary_path = (
                            sec_path_uv if sec_u == hi_atom else list(reversed(sec_path_uv))
                        )
                        secondary.append(
                            {"lo": lo, "hi": hi, "path": secondary_path, "len": len(sec_path_uv)}
                        )
                    if not all_independent:
                        continue

                    # VB-7: independent secondary bridges are numbered in
                    # order of the highest-numbered bridgehead they attach
                    # to (higher first), ties broken by longer bridge, then
                    # lower attachment locant.
                    numbering_order = sorted(
                        secondary, key=lambda s: (-s["hi"], -s["len"], s["lo"])
                    )
                    full_order = list(main_order)
                    for s in numbering_order:
                        full_order.extend(s["path"])

                    # VB-6: secondary bridges are cited in decreasing size;
                    # equal-size ties cited in ascending locant order (this is
                    # independent of, and can differ from, numbering_order --
                    # see the tetracyclo[4.4.2.2^2,5.2^7,10]hexadecane example
                    # in this module's docstring).
                    citation_order = sorted(
                        secondary, key=lambda s: (-s["len"], s["lo"], s["hi"])
                    )
                    descriptor = ".".join(
                        f"{s['len']}^{s['lo']},{s['hi']}" for s in citation_order
                    )
                    total_atoms = a + b + c + sum(s["len"] for s in secondary) + 2
                    parent = f"{ring_prefix}[{a}.{b}.{c}.{descriptor}]{alkane_name(total_atoms)}"

                    combined_locants = tuple(
                        sorted(loc for s in secondary for loc in (s["lo"], s["hi"]))
                    )
                    citation_locants = tuple(
                        loc for s in citation_order for loc in (s["lo"], s["hi"])
                    )
                    # VB-6 (maximize each secondary bridge's length, largest
                    # first) > P-23.2.6.2.4/.2.5-style lowest combined locant
                    # set for all secondary-bridge attachment points > lowest
                    # locants in citation order > (caller's own tie-break,
                    # e.g. P-14.4/P-45.2 lowest substituent locants).
                    outer_key = (
                        shape_key
                        + tuple(-s["len"] for s in citation_order)
                        + combined_locants
                        + citation_locants
                    )
                    yield full_order, parent, outer_key


def _candidate_key(parent, substituents, heteroatom_locant=None, nondetachable_prefix=""):
    """`heteroatom_locant`/`nondetachable_prefix`: as of
    `tasks/heteroatom-skeleton-expansion.md` (2026-08-26), shared with
    `_von_baeyer_heteroatom.py`'s polycyclic (ring_count>=3) case, mirroring
    `_bicyclic._candidate_key`'s same two optional parameters -- see that
    function's callers for why the heteroatom locant is ranked ahead of
    substituent locants but after `iter_polycyclic_candidates`'s own
    `outer_key` (which the caller combines with this function's return
    value)."""
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    if prefix and nondetachable_prefix:
        prefix += "-"
    name = prefix + nondetachable_prefix + parent
    if heteroatom_locant is None:
        return locant_set, citation_locants, name
    return heteroatom_locant, locant_set, citation_locants, name


def name_polycycloalkane(mol, core, ring_count) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated polycyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    best_name = None
    for full_order, parent, outer_key in iter_polycyclic_candidates(core, ring_count):
        substituents = _substituents_for_ring(graph, full_order, halogens)
        key = outer_key + _candidate_key(parent, substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, key[-1]

    if best_name is None:
        raise UnsupportedStructure(
            "this polycyclic topology is not supported yet (see "
            "_polycyclic.py's module docstring for the scope this module "
            "covers)"
        )
    return best_name
