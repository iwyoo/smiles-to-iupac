"""Von Baeyer naming of saturated tetracyclic hydrocarbons whose skeleton
reduces to exactly six skeletal atoms of degree 3 (cyclomatic number 4: a
main bicyclic system per P-23.2.2/P-23.2.3 -- two main bridgeheads and three
bridges between them -- plus two further *independent* secondary bridges
among the remaining four branch atoms, each attaching only to atoms already
placed by the main system, per the same "independent secondary bridge"
notion P-23.2.5.1 uses for the tricyclic case, just doubled).

This generalizes `_tricyclic.py`'s "pick 3 composite bridges between 2 main
bridgeheads, 1 bridge left over as the secondary bridge" search by one more
independent bridge: cyclomatic number 3 -> 4 means one more edge than
vertices-implied-by-a-tree, i.e. one more independent bridge. Six branch
atoms of degree 3 give 6*3/2 = 9 real bridges between them (vs. tricyclic's
four branch atoms / six bridges); three of those nine form the main system
(as in `_bicyclic.py`/`_tricyclic.py`) and the other six real bridges must
resolve into exactly two independent secondary bridges plus four "hop"
edges that place the four non-main-bridgehead branch atoms onto the main
system's own numbering (each such branch atom uses two of its three edges
as an intermediate hop on one of the three main-system composite bridges,
and its third edge as one endpoint of one of the two secondary bridges --
see `_composite_bridges` below).

Propellane-like topologies (any branch atom of degree 4 or higher) and
pentacyclic-or-higher systems (cyclomatic number >= 5) remain out of scope
and raise UnsupportedStructure, exactly as degree-4 atoms and
tetracyclic-or-higher systems do in `_tricyclic.py` -- see
`find_tetracyclic_core`'s degree/branch-atom-count checks. So does any
tetracyclic topology that doesn't reduce to six branch atoms with two
*independent* secondary bridges (e.g. a "dependent" secondary bridge whose
attachment point lies on another secondary bridge rather than on the main
ring/main bridge itself) -- `name_tetracycloalkane` simply finds no valid
decomposition for such cases and raises.

Per the IUPAC 2013 Recommendations ("the Blue Book", Chapter P-2,
https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf), P-23.2.6 ("polycyclic ring
systems containing four or more rings") extends P-23.2.1/P-23.2.4/P-23.2.5
by incorporating, largely unchanged, the earlier IUPAC "von Baeyer" rules
VB-1 through VB-9 (W.H. Powell, "Extension and Revision of the von Baeyer
System for Naming Polycyclic Compounds", Pure Appl. Chem. 71 (1999)
513-529; the numbered VB rules are mirrored at
https://iupac.qmul.ac.uk/vonBaeyer/vb6.html and .../vb7.html, which is what
this module's multi-secondary-bridge handling was verified against, since
the Blue Book's own P-23.2.6.x text was not text-extractable from its PDF
at the time this module was written -- see this module's final
implementation report for that caveat):

- P-23.2.1 (max main ring) and P-23.2.4 (max main bridge): unchanged from
  the bicyclic/tricyclic case (see `_bicyclic.py`/`_tricyclic.py`); still
  take priority over everything below, per P-23.2.6.2.1's incorporation of
  those same two rules for the tetracyclic-and-higher case.
- VB-6 ("Naming Hydrocarbon Polycyclic Systems"): "The numbers indicating
  independent secondary bridges are cited in decreasing order" -- e.g. the
  worked example there prefers tetracyclo[8.6.6.5^2,9.1^23,26]octacosane
  over tetracyclo[8.6.6.4^2,9.2^23,25]octacosane ("8.6.6.5.1 is higher than
  8.6.6.4.2"): given a choice of which real bridges become the two
  secondary bridges (after the main ring and main bridge are fixed), the
  decomposition maximizing the *larger* secondary bridge's length wins,
  then the smaller one -- applied here as the next tie-break after
  P-23.2.1/P-23.2.4, mirroring how P-23.2.4 itself maximizes the main
  bridge only after the main ring is already fixed.
- VB-7 ("Numbering of Secondary Bridges"): "After numbering the main ring
  and main bridge, all independent secondary bridges are numbered before
  dependent secondary bridges" (not reached here -- this module only
  supports independent/independent pairs), continuing "from the highest
  number of the main ring and bridge", with the secondary bridge *attached
  to the highest-numbered bridgehead numbered first* -- i.e. whichever of
  the two secondary bridges' higher attachment locant is larger has its own
  internal atoms numbered first (gets the lower of the two "new" locant
  ranges). VB-7.1: ties broken by giving lower locants to the bridge linked
  to the higher-numbered bridgehead; VB-7.2: ties still remaining broken by
  numbering the longer bridge first. (Within a single secondary bridge,
  P-23.2.5.1/.2's already-established convention is kept: numbering starts
  adjacent to the *higher*-numbered attachment bridgehead, but the pair is
  still *cited* lower-locant-first.)
- Citation order in the descriptor is a separate concern from numbering
  order above: secondary bridges are cited in decreasing size (VB-6, see
  above); the tetracyclo[4.4.2.2^2,5.2^7,10]hexadecane worked example
  (VB-7B) has two equal-length (2-atom) secondary bridges, cited in
  ascending order of their locant pairs (2,5 before 7,10) even though the
  (7,10) bridge is *numbered* first per VB-7 -- confirming citation order
  and numbering order are independently determined.
- P-23.2.6.2.4/.2.5-style tie-breaks (mirroring `_tricyclic.py`'s own
  generalization of these for its single secondary bridge): once main
  ring/main bridge and the secondary bridges' own lengths are fixed, prefer
  the decomposition with the lowest combined locant set for all four
  secondary-bridge attachment points, then the lowest locants in citation
  order.
- P-14.4 / P-45.2 (Chapter P-1 / P-4): any remaining choice is settled by
  lowest locants to substituents, exactly as `_bicyclic.py`/`_tricyclic.py`
  do.

Flagship validation case: quadricyclane, tetracyclo[3.2.0.0^2,7.0^4,6]heptane
(C7H8; PubChem CID 78961 gives canonical SMILES C1C2C3C2C4C1C34 and this
exact PubChem-computed IUPAC name -- cross-checked externally via the
PubChem PUG REST API, not from memory, and independently re-verified via
RDKit's reported formula/degree sequence: six atoms of degree 3, one of
degree 2, consistent with this module's six-branch-atom scope).

Fused/spiro/bicyclic/tricyclic systems are handled by
`_cyclic.py`/`_spiro.py`/`_bicyclic.py`/`_tricyclic.py`; unsaturated and
heteroatom-containing tetracyclics, pentacyclic-or-higher systems, and any
tetracyclic topology this module's detection can't cleanly resolve into the
scope above (propellane-like degree-4 branch atoms, or dependent secondary
bridges) raise UnsupportedStructure.
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


def find_tetracyclic_core(mol):
    """Return (branch_atoms, bridges) if `mol`'s carbon skeleton, after
    stripping acyclic branches, reduces to exactly six branch atoms of
    degree 3 (and no atom of higher degree) joined by nine bridges (the
    graph's cyclomatic number is 4) -- else None.

    `bridges` is a list of nine `(u, v, path)` tuples, one per bridge, where
    `path` is the list of internal atom indices ordered nearest-`u`-first.
    As in `_tricyclic.find_tricyclic_core`, branch-atom pairs need not each
    have exactly one bridge between them -- doubled bridges between the
    same pair are a real, connected topology and are detected the same way.

    Mirrors `_tricyclic.find_tricyclic_core`'s leaf-stripping +
    cyclomatic-number approach (see that function's docstring for why this
    is preferred over trusting RDKit's `NumRings()`)."""
    graph = adjacency(mol)
    core = _strip_leaves(graph)
    if not core:
        return None
    vertices = len(core)
    edge_count = sum(len(neighbors) for neighbors in core.values()) // 2
    if edge_count - vertices + 1 != 4:
        return None
    if any(len(neighbors) not in (2, 3) for neighbors in core.values()):
        return None
    branch_atoms = {atom for atom, neighbors in core.items() if len(neighbors) == 3}
    if len(branch_atoms) != 6:
        return None

    half_walks = []
    for u in branch_atoms:
        for first in core[u]:
            result = _walk_to_branch(core, branch_atoms, u, first)
            if result is None:
                return None
            v, path = result
            half_walks.append((u, v, path))
    if len(half_walks) != 18:
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
    if len(bridges) != 9:
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
    any intermediate branch atom(s).

    Generalizes `_tricyclic._composite_bridges` (which only ever has to
    consider zero, one, or both of exactly *two* "other" branch atoms as
    intermediate hops, since a tricyclic system has only four branch atoms
    total) to an arbitrary number of "other" atoms -- the tetracyclic case
    has four remaining branch atoms rather than two, so a single composite
    bridge here may need to hop through up to all four of them (see
    quadricyclane's main-ring-second segment in this module's flagship test,
    which hops through two of the four). `_tricyclic._composite_bridges`'s
    fixed signature (`others` unpacked as exactly `o1, o2 = others`) doesn't
    generalize to a variable-size `others`, so this is a standalone
    reimplementation of the same idea (a bounded-depth backtracking walk)
    rather than a direct import."""
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


def name_tetracycloalkane(mol, core) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated tetracyclic ring systems are not supported yet (see "
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
                if len(remaining) != 2:
                    continue
                sec_idx1, sec_idx2 = sorted(remaining)

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

                    secondary = []
                    both_independent = True
                    for sec_idx in (sec_idx1, sec_idx2):
                        sec_u, sec_v, sec_path_uv = bridges[sec_idx]
                        if sec_u not in position or sec_v not in position:
                            both_independent = False
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
                    if not both_independent:
                        continue

                    # VB-7: the independent secondary bridge attached to the
                    # highest-numbered bridgehead is numbered first (its own
                    # atoms get the lower of the two new locant ranges);
                    # VB-7.2/VB-7.1-style ties broken by longer bridge, then
                    # lower attachment locant.
                    numbering_order = sorted(secondary, key=lambda s: (-s["hi"], -s["len"], s["lo"]))
                    full_order = list(main_order)
                    for s in numbering_order:
                        full_order.extend(s["path"])

                    substituents = _substituents_for_ring(graph, full_order, halogens)

                    # VB-6: secondary bridges are cited in decreasing size;
                    # equal-size ties cited in ascending locant order (this is
                    # independent of, and can differ from, numbering_order --
                    # see the tetracyclo[4.4.2.2^2,5.2^7,10]hexadecane example
                    # in this module's docstring).
                    citation_order = sorted(secondary, key=lambda s: (-s["len"], s["lo"], s["hi"]))
                    d1, lo1, hi1 = (
                        citation_order[0]["len"],
                        citation_order[0]["lo"],
                        citation_order[0]["hi"],
                    )
                    d2, lo2, hi2 = (
                        citation_order[1]["len"],
                        citation_order[1]["lo"],
                        citation_order[1]["hi"],
                    )

                    total_atoms = a + b + c + d1 + d2 + 2
                    parent = (
                        f"tetracyclo[{a}.{b}.{c}.{d1}^{lo1},{hi1}.{d2}^{lo2},{hi2}]"
                        f"{alkane_name(total_atoms)}"
                    )
                    # VB-6 (maximize the larger secondary bridge, then the
                    # smaller) > P-23.2.6.2.4/.2.5-style lowest combined
                    # locant set for all secondary-bridge attachment points >
                    # lowest locants in citation order > P-14.4/P-45.2 lowest
                    # substituent locants.
                    key = (
                        outer_key
                        + (-d1, -d2)
                        + tuple(sorted((lo1, hi1, lo2, hi2)))
                        + (lo1, hi1, lo2, hi2)
                        + _candidate_key(parent, substituents)
                    )
                    if best_key is None or key < best_key:
                        best_key, best_name = key, key[-1]

    if best_name is None:
        raise UnsupportedStructure(
            "this tetracyclic topology is not supported yet (see "
            "_tetracyclic.py's module docstring for the scope this module "
            "covers)"
        )
    return best_name
