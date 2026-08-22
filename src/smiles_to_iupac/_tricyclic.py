"""Von Baeyer naming of saturated tricyclic hydrocarbons whose three rings
reduce to a bicyclic main system (two main bridgeheads, three bridges, as in
`_bicyclic.py`) plus one independent secondary bridge connecting two further,
distinct atoms that already lie on that main bicyclic system (the "secondary
bridgeheads") -- i.e. exactly four skeletal atoms of degree 3 and none of
higher degree. Other tricyclic topologies (e.g. a secondary bridge that
reconnects to a main bridgehead itself, collapsing the four branch points to
fewer distinct atoms -- a propellane-like case) and tetracyclic-or-higher
systems are out of scope and raise UnsupportedStructure. Per the IUPAC 2013
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
  main bridge, exactly as `_bicyclic.py` computes them; d the secondary
  bridge's length; x,y its attachment locants, lower first, rendered here as
  plain-text `d^x,y` with no braces) + the alkane name for the total number
  of skeletal atoms.
- P-23.2.5.2: "After the main ring and main bridge have been numbered, the
  independent secondary bridge is numbered continuing from the higher
  numbered bridgehead of the main ring." (confirmed by the worked example
  "tricyclo[4.2.2.2^2,5]dodecane [the secondary bridge is numbered starting
  from bridgehead 5, the higher (than 2) numbered bridgehead]" -- numbering
  starts adjacent to the HIGHER-locant secondary bridgehead, not the lower
  one, even though the lower one is cited first in the bracket descriptor.)
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

Flagship validation case: adamantane, tricyclo[3.3.1.1^3,7]decane (verified
structurally here -- C10H16, four skeletal atoms of degree 3 and six of
degree 2 -- and cross-checked against the NIST WebBook / PubChem name for
CAS 281-23-2, and against the derivation above applied by hand).

Fused/spiro/bicyclic systems are handled by `_cyclic.py`/`_spiro.py`/
`_bicyclic.py`; unsaturated and heteroatom-containing tricyclics, and any
tricyclic topology this module's detection can't cleanly resolve into the
scope above, raise UnsupportedStructure.
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


def find_tricyclic_core(mol):
    """Return (branch_atoms, edges) if `mol`'s carbon skeleton, after
    stripping acyclic branches, reduces to exactly four branch atoms of
    degree 3 (and no atom of higher degree), pairwise connected by six
    bridges -- i.e. a complete graph K4 on the four branch atoms, the
    topology of a bicyclic main system plus one independent secondary bridge
    between two of its non-main-bridgehead atoms (see module docstring) --
    else None.

    `edges` is {frozenset({u, v}): {u: path_from_u, v: path_from_v}} for each
    of the six branch-atom pairs; each path is the list of internal atom
    indices ordered starting next to the named endpoint.

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

    edges = {}
    for u in branch_atoms:
        for first in core[u]:
            result = _walk_to_branch(core, branch_atoms, u, first)
            if result is None:
                return None
            v, path = result
            sides = edges.setdefault(frozenset((u, v)), {})
            if u in sides:
                return None
            sides[u] = path
    if len(edges) != 6:
        return None
    for sides in edges.values():
        if len(sides) != 2:
            return None
        (u, path_u), (v, path_v) = sides.items()
        if path_u != list(reversed(path_v)):
            return None

    return branch_atoms, edges


def _bridge_from(edges, u, v):
    return edges[frozenset((u, v))][u]


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

    branch_atoms, edges = core
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    best_name = None
    for bh1, bh2 in combinations(branch_atoms, 2):
        sec1, sec2 = (a for a in branch_atoms if a not in (bh1, bh2))
        bridge_direct = _bridge_from(edges, bh1, bh2)
        bridge_via_sec1 = _bridge_from(edges, bh1, sec1) + [sec1] + _bridge_from(edges, sec1, bh2)
        bridge_via_sec2 = _bridge_from(edges, bh1, sec2) + [sec2] + _bridge_from(edges, sec2, bh2)
        bridges = [bridge_direct, bridge_via_sec1, bridge_via_sec2]

        a, b, c = sorted((len(bridge) for bridge in bridges), reverse=True)
        d = len(_bridge_from(edges, sec1, sec2))
        total_atoms = a + b + c + d + 2
        # P-23.2.1 (max main ring) > P-23.2.4 (max main bridge) >
        # P-23.2.6.2.1 (main ring divided as symmetrically as possible).
        outer_key = (-(a + b), -c, a - b)

        for start, other, oriented_bridges in (
            (bh1, bh2, bridges),
            (bh2, bh1, [list(reversed(bridge)) for bridge in bridges]),
        ):
            for perm in permutations(range(3)):
                ordered = [oriented_bridges[i] for i in perm]
                if not (len(ordered[0]) >= len(ordered[1]) >= len(ordered[2])):
                    continue
                main_ring_first, main_ring_second, main_bridge = ordered
                main_order = (
                    [start] + main_ring_first + [other] + list(reversed(main_ring_second)) + main_bridge
                )
                position = {atom: i + 1 for i, atom in enumerate(main_order)}
                x, y = position[sec1], position[sec2]
                lo, hi = min(x, y), max(x, y)
                hi_atom, lo_atom = (sec1, sec2) if x > y else (sec2, sec1)
                secondary_path = _bridge_from(edges, hi_atom, lo_atom)
                full_order = main_order + secondary_path

                substituents = _substituents_for_ring(graph, full_order, halogens)
                parent = f"tricyclo[{a}.{b}.{c}.{d}^{lo},{hi}]{alkane_name(total_atoms)}"
                # P-23.2.6.2.4/.2.5 (lowest secondary-bridge locants) then
                # P-14.4/P-45.2 (lowest substituent locants).
                key = outer_key + ((lo, hi),) + _candidate_key(parent, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]

    return best_name
