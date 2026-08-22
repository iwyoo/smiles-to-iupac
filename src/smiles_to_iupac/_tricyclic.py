"""Von Baeyer naming of saturated tricyclic hydrocarbons whose skeleton
reduces to exactly four skeletal atoms of degree 3 (a "main ring" plus a
main bridge plus one independent secondary bridge, per P-23.2.5.1) and none
of higher degree. This covers two structurally distinct branch-atom
multigraphs:

- the "K4" case, where each of the four branch atoms connects directly (via
  a single bridge with no other branch atom on it) to each of the other
  three -- e.g. adamantane, twistane;
- the "doubled main bridgeheads" case, where two branch atoms are joined by
  *two* parallel bridges and the other two branch atoms are also joined by
  two parallel bridges, with a single bridge connecting one atom from each
  pair -- this is the Blue Book's own P-23.2.5.2 worked example,
  tricyclo[4.2.2.2^2,5]dodecane, where the two secondary-bridge attachment
  points (locants 2 and 5) fall on the same 4-atom main-ring segment. It
  also covers ortho-fused "chain of rings" tricyclics such as
  perhydroanthracene/perhydrophenanthrene, which are the same abstract
  branch-atom multigraph with different bridge lengths.

Propellane-like topologies -- two branch atoms, both of degree 4, directly
bonded to each other *in addition to* three bridges between them -- are a
structurally distinct case, handled separately by `find_propellane_core`/
`name_propellane` below rather than by the four-branch-atom search above:
the direct bond is itself the independent secondary bridge of P-23.2.5.1,
with length 0 (`0^x,y`, same zero-length-bridge notation already used
elsewhere in this module, e.g. tricyclo[4.4.0.0^3,8]decane) and its two
attachment points are the main bridgeheads themselves. Tetracyclic-or-higher
systems remain out of scope and raise UnsupportedStructure -- see
`find_tricyclic_core`'s and `find_propellane_core`'s degree/branch-atom-count
checks. Per the IUPAC 2013
Recommendations ("the Blue Book", Chapter P-2,
https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):

- P-23.2.1: "The main ring of a polycyclic hydrocarbon ring system is
  selected so as to include as many skeletal atoms of the structure as
  possible."
- P-23.2.4: "...the main bridge...is the bridge that includes as many of the
  atoms as possible that are not included in the main ring. Bridges other
  than the main bridge are called 'secondary bridges'."
- P-23.2.5.1: "Tricyclic hydrocarbons having an independent secondary bridge
  are named on the basis of a bicyclic system described in P-23.2.2. Rings
  not described by the bicyclic system are defined by citing the number of
  atoms in the independent secondary bridge as an arabic number. The locants
  of the two attachment points of the independent secondary bridge to the
  main ring are cited as a pair of superscript arabic numbers (lower number
  is cited first) separated by a comma." The name is 'tricyclo' + bridge
  lengths `[a.b.c.d^x,y]` (a >= b >= c the two main-ring segments and the
  main bridge; d the secondary bridge's length; x,y its attachment locants,
  lower first, rendered here as plain-text `d^x,y` with no braces) + the
  alkane name for the total number of skeletal atoms.
- P-23.2.5.2: "After the main ring and main bridge have been numbered, the
  independent secondary bridge is numbered continuing from the higher
  numbered bridgehead of the main ring." (confirmed by the worked example
  "tricyclo[4.2.2.2^2,5]dodecane [the secondary bridge is numbered starting
  from bridgehead 5, the higher (than 2) numbered bridgehead]" -- numbering
  starts adjacent to the HIGHER-locant secondary bridgehead, not the lower
  one, even though the lower one is cited first in the bracket descriptor.
  This worked example is itself the "doubled main bridgeheads" case above:
  both attachment points 2 and 5 lie on the same 4-atom main-ring segment,
  meaning the main bridgeheads are joined by two parallel bridges (the two
  length-2 segments b, c) and the secondary bridgeheads -- interior points
  of the length-4 segment a -- are joined by two parallel bridges as well
  (the interior of a itself, and the independent secondary bridge d).)
- P-23.2.6.2.1/.2.4/.2.5 (stated for tetracyclic-and-higher systems, but
  generalizing P-23.2.1/P-23.2.4's selection criteria to any case where more
  than one choice of main bridgeheads is possible, so applied here too):
  after maximizing the main ring (P-23.2.1) and then the main bridge
  (P-23.2.4), remaining ties are broken by (a) dividing the main ring as
  symmetrically as possible between its two segments, then (b) the lowest
  possible locant set for the secondary bridge's attachment points, then (c)
  the lowest locants in citation order.
- P-14.4 / P-45.2 (Chapter P-1 / P-4): any choice the rules above still leave
  open (which main bridgehead is numbered first, which of two equal-length
  bridges is traversed first) is settled by lowest locants to substituents,
  exactly as `_bicyclic.py` does.

Flagship validation cases:
- adamantane, tricyclo[3.3.1.1^3,7]decane (K4 case; C10H16, four skeletal
  atoms of degree 3 and six of degree 2 -- cross-checked against the NIST
  WebBook / PubChem name for CAS 281-23-2).
- tricyclo[4.2.2.2^2,5]dodecane (doubled-main-bridgeheads case; built
  directly from the Blue Book's own P-23.2.5.2 worked example, C12H20).
- [1.1.1]propellane, tricyclo[1.1.1.0^1,3]pentane (propellane case; C5H6 --
  cross-checked against Wikipedia/ChemSpider/ACS "Molecule of the Week").
- [2.2.2]propellane, tricyclo[2.2.2.0^1,4]octane (propellane case; C8H12 --
  cross-checked against Wikipedia/Wikidata).

Fused/spiro/bicyclic systems are handled by `_cyclic.py`/`_spiro.py`/
`_bicyclic.py`; unsaturated and heteroatom-containing tricyclics, and any
tricyclic topology this module's detection can't cleanly resolve into the
scope above (propellanes, tetracyclic-or-higher), raise UnsupportedStructure.
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
from ._numerals import alkane_name
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


def find_tricyclic_core(mol):
    """Return (branch_atoms, bridges) if `mol`'s carbon skeleton, after
    stripping acyclic branches, reduces to exactly four branch atoms of
    degree 3 (and no atom of higher degree) joined by six bridges (the
    graph's cyclomatic number is 3) -- else None.

    `bridges` is a list of six `(u, v, path)` tuples, one per bridge, where
    `path` is the list of internal atom indices ordered nearest-`u`-first.
    Unlike the topology this module originally supported, branch-atom pairs
    need not each have exactly one bridge between them: two bridges between
    the same pair (and, correspondingly, some other pair having none) is a
    real, connected topology -- see the module docstring -- and is detected
    here the same way.

    Mirrors `_bicyclic.find_bicyclic_core`'s leaf-stripping + cyclomatic-
    number approach rather than trusting RDKit's `NumRings()` (see that
    function's docstring for why)."""
    graph = adjacency(mol)
    core = _strip_leaves(graph)
    if not core:
        return None
    vertices = len(core)
    edge_count = sum(len(neighbors) for neighbors in core.values()) // 2
    if edge_count - vertices + 1 != 3:
        return None
    if any(len(neighbors) not in (2, 3) for neighbors in core.values()):
        return None
    branch_atoms = {atom for atom, neighbors in core.items() if len(neighbors) == 3}
    if len(branch_atoms) != 4:
        return None

    half_walks = []
    for u in branch_atoms:
        for first in core[u]:
            result = _walk_to_branch(core, branch_atoms, u, first)
            if result is None:
                return None
            v, path = result
            half_walks.append((u, v, path))
    if len(half_walks) != 12:
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
    if len(bridges) != 6:
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
    `start` to `end` that uses zero, one, or both of `others` as
    intermediate branch atoms, each hop a real bridge with no bridge index
    reused. `atom_path` is the full ordered list of atoms strictly between
    `start` and `end`, including any intermediate branch atom(s) -- so when
    a composite is later used as a main-system bridge, its intermediate
    branch atom(s) fall naturally into the resulting numbering."""
    o1, o2 = others
    for idx, v, path in adj[start]:
        if v == end:
            yield (idx,), path
    for mid in (o1, o2):
        for idx1, v1, path1 in adj[start]:
            if v1 != mid:
                continue
            for idx2, v2, path2 in adj[mid]:
                if v2 != end or idx2 == idx1:
                    continue
                yield (idx1, idx2), path1 + [mid] + path2
    for mid1, mid2 in ((o1, o2), (o2, o1)):
        for idx1, v1, path1 in adj[start]:
            if v1 != mid1:
                continue
            for idx2, v2, path2 in adj[mid1]:
                if v2 != mid2 or idx2 == idx1:
                    continue
                for idx3, v3, path3 in adj[mid2]:
                    if v3 != end or idx3 in (idx1, idx2):
                        continue
                    yield (
                        (idx1, idx2, idx3),
                        path1 + [mid1] + path2 + [mid2] + path3,
                    )


def _candidate_key(parent, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = parent if not grouped else format_substituent_prefixes(grouped) + parent
    return locant_set, citation_locants, name


def name_tricycloalkane(mol, core) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated tricyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    branch_atoms, bridges = core
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    adj = _adjacency_by_atom(branch_atoms, bridges)

    best_key = None
    best_name = None
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
                if len(remaining) != 1:
                    continue
                sec_idx = next(iter(remaining))
                sec_u, sec_v, sec_path_uv = bridges[sec_idx]

                a, b, c = sorted((len(p) for p in paths), reverse=True)
                # P-23.2.1 (max main ring) > P-23.2.4 (max main bridge) >
                # P-23.2.6.2.1 (main ring divided as symmetrically as possible).
                outer_key = (-(a + b), -c, a - b)

                for perm in permutations(paths):
                    if not (len(perm[0]) >= len(perm[1]) >= len(perm[2])):
                        continue
                    main_ring_first, main_ring_second, main_bridge = perm
                    main_order = (
                        [start] + main_ring_first + [other]
                        + list(reversed(main_ring_second)) + main_bridge
                    )
                    position = {atom: idx + 1 for idx, atom in enumerate(main_order)}
                    if sec_u not in position or sec_v not in position:
                        continue
                    x, y = position[sec_u], position[sec_v]
                    lo, hi = min(x, y), max(x, y)
                    hi_atom = sec_u if x > y else sec_v
                    secondary_path = (
                        sec_path_uv if sec_u == hi_atom else list(reversed(sec_path_uv))
                    )
                    full_order = main_order + secondary_path

                    substituents = _substituents_for_ring(graph, full_order, halogens)
                    d = len(sec_path_uv)
                    total_atoms = a + b + c + d + 2
                    parent = f"tricyclo[{a}.{b}.{c}.{d}^{lo},{hi}]{alkane_name(total_atoms)}"
                    # P-23.2.6.2.4/.2.5 (lowest secondary-bridge locants) then
                    # P-14.4/P-45.2 (lowest substituent locants).
                    key = outer_key + ((lo, hi),) + _candidate_key(parent, substituents)
                    if best_key is None or key < best_key:
                        best_key, best_name = key, key[-1]

    if best_name is None:
        raise UnsupportedStructure(
            "this tricyclic topology is not supported yet (see "
            "_tricyclic.py's module docstring for the scope this module "
            "covers)"
        )
    return best_name


def find_propellane_core(mol):
    """Return (bh1, bh2, bridges) if `mol`'s carbon skeleton, after stripping
    acyclic branches, reduces to exactly two branch atoms of degree 4 (and no
    atom of higher degree), directly bonded to each other, joined by exactly
    three further bridges through degree-2 atoms -- else None.

    `bridges` is a list of three `(bh1, bh2, path)` tuples, `path` ordered
    nearest-`bh1`-first. The cyclomatic-number-3 check makes this genuinely
    tricyclic (not tetracyclic-or-higher) for any bridge lengths: with two
    degree-4 branch atoms and three bridges plus the direct bond, edges -
    vertices + 1 is always 3 regardless of bridge length, so no extra guard
    against tetracyclic-or-higher is needed beyond the checks already here."""
    graph = adjacency(mol)
    core = _strip_leaves(graph)
    if not core:
        return None
    vertices = len(core)
    edge_count = sum(len(neighbors) for neighbors in core.values()) // 2
    if edge_count - vertices + 1 != 3:
        return None
    if any(len(neighbors) not in (2, 4) for neighbors in core.values()):
        return None
    branch_atoms = {atom for atom, neighbors in core.items() if len(neighbors) == 4}
    if len(branch_atoms) != 2:
        return None
    bh1, bh2 = branch_atoms
    if bh2 not in core[bh1]:
        return None

    bridges = []
    for first in core[bh1] - {bh2}:
        result = _walk_to_branch(core, branch_atoms, bh1, first)
        if result is None:
            return None
        v, path = result
        if v != bh2:
            return None
        bridges.append((bh1, bh2, path))
    if len(bridges) != 3:
        return None

    return bh1, bh2, bridges


def name_propellane(mol, core) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated propellane-type ring systems are not supported yet "
            "(see P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    bh1, bh2, bridges = core
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    best_name = None
    for start, other in ((bh1, bh2), (bh2, bh1)):
        oriented_paths = [
            path if start == bh1 else list(reversed(path)) for _, _, path in bridges
        ]
        for perm in permutations(oriented_paths):
            if not (len(perm[0]) >= len(perm[1]) >= len(perm[2])):
                continue
            main_ring_first, main_ring_second, main_bridge = perm
            main_order = (
                [start] + list(main_ring_first) + [other]
                + list(reversed(main_ring_second)) + list(main_bridge)
            )
            position = {atom: idx + 1 for idx, atom in enumerate(main_order)}
            lo, hi = sorted((position[start], position[other]))
            substituents = _substituents_for_ring(graph, main_order, halogens)
            a, b, c = len(main_ring_first), len(main_ring_second), len(main_bridge)
            total_atoms = a + b + c + 2
            # The direct bond between the two bridgeheads is the independent
            # secondary bridge (P-23.2.5.1), length 0, attached at the main
            # bridgeheads themselves -- see module docstring.
            parent = f"tricyclo[{a}.{b}.{c}.0^{lo},{hi}]{alkane_name(total_atoms)}"
            # P-23.2.6.2.4/.2.5 (lowest secondary-bridge locants) then
            # P-14.4/P-45.2 (lowest substituent locants).
            key = ((lo, hi),) + _candidate_key(parent, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    if best_name is None:
        raise UnsupportedStructure(
            "this propellane topology is not supported yet (see "
            "_tricyclic.py's module docstring for the scope this module "
            "covers)"
        )
    return best_name
