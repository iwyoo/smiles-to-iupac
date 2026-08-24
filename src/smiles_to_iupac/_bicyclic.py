"""Von Baeyer naming of saturated bicyclic hydrocarbons (two carbocyclic rings
sharing two or more atoms, fused or bridged), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-23.1.1 / P-23.1.2 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  a 'bridgehead' is any skeletal ring atom bonded to three or more skeletal ring
  atoms; a 'bridge' is an unbranched chain of atoms, a single atom, or a valence
  bond, connecting two bridgeheads. The systems handled here have exactly two
  bridgeheads and three bridges between them (one bridge may have zero atoms,
  i.e. the bridgeheads directly bonded, as in a fused ring system).
- P-23.2.2: "Saturated homogeneous bicyclic hydrocarbons having two or more
  atoms in common are named by prefixing 'bicyclo' to the name of the acyclic
  hydrocarbon having the same total number of skeletal atoms. The numbers of
  skeletal atoms in each of the two segments connecting the main bridgeheads
  and in the main bridge are given by arabic numbers cited in descending
  numerical order, separated by full stops, and enclosed in square brackets."
  (e.g. bicyclo[3.2.1]octane, bicyclo[2.2.1]heptane, bicyclo[4.4.0]decane).
- P-23.2.3: "The bicyclic ring system is numbered starting with one of the
  bridgeheads and proceeding first along the longer segment of the main ring
  to the second bridgehead, then back to the first bridgehead along the
  unnumbered segment of the main ring. Numbering is completed by numbering the
  main bridge beginning with the atom next to the first bridgehead."
- P-14.4 / P-45.2 (Chapter P-1 / P-4): when P-23.2.3 leaves a choice (which
  bridgehead is numbered first, which of two equal-length bridges is
  numbered before the other), lowest locants go to substituents as a set,
  then in order of citation — the same tie-break machinery `_cyclic.py` and
  `_spiro.py` use.
- P-29.4 / P-46: branched ("compound") substituent groups — see
  `_substituents.py`.

- P-23.0: "This Section deals only with saturated polyalicyclic ring systems
  named by the von Baeyer system; for unsaturated systems, see Section
  P-31.1.4."

Tricyclic and higher polycyclic systems (P-23.2.5, P-23.2.6), unsaturated
bicyclics (P-31.1.4), and heteroatom-containing bicyclics (P-23.3) are out of
scope and raise UnsupportedStructure.
"""

from itertools import permutations

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
    ring-system core (see module docstring's bridgehead/bridge definitions)."""
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


def _walk_bridge(core, bh1, bh2, first):
    """Walk the bridge starting at `first` (a neighbor of bh1) toward bh2,
    returning its internal atoms in order (nearest-bh1 first), or None if the
    walk doesn't cleanly reach bh2 (e.g. two separate rings joined by a single
    bond, which share no atom and so aren't a bicyclic system at all)."""
    if first == bh2:
        return []
    path = [first]
    previous, current = bh1, first
    while True:
        neighbors = core[current] - {previous}
        if len(neighbors) != 1:
            return None
        next_atom = next(iter(neighbors))
        if next_atom == bh2:
            return path
        if next_atom == bh1:
            return None
        path.append(next_atom)
        previous, current = current, next_atom


def find_bicyclic_core(mol):
    """Return (bridgehead1, bridgehead2, [bridge1, bridge2, bridge3]) if `mol`'s
    carbon skeleton, after stripping acyclic branches, reduces to exactly two
    bridgeheads joined by three bridges (P-23.1.1/P-23.1.2's scope: a
    genuinely bicyclic hydrocarbon), else None. Each bridge is the list of its
    internal atom indices, ordered from the bh1 side to the bh2 side.

    Deliberately does not rely on RDKit's `GetRingInfo().NumRings()`: its SSSR
    can return a redundant, non-minimal ring set for symmetric bridged
    bicyclics (e.g. bicyclo[2.2.2]octane reports 3 rings, not 2), so ring
    count is instead derived from the core's own cyclomatic number
    (edges - vertices + 1), which is exactly 2 for a genuine bicyclic core."""
    graph = adjacency(mol)
    core = _strip_leaves(graph)
    if not core:
        return None
    vertices = len(core)
    edges = sum(len(neighbors) for neighbors in core.values()) // 2
    if edges - vertices + 1 != 2:
        return None
    if any(len(neighbors) not in (2, 3) for neighbors in core.values()):
        return None
    bridgeheads = [atom for atom, neighbors in core.items() if len(neighbors) == 3]
    if len(bridgeheads) != 2:
        return None
    bh1, bh2 = bridgeheads
    bridges = [_walk_bridge(core, bh1, bh2, start) for start in core[bh1]]
    if any(bridge is None for bridge in bridges):
        return None
    return bh1, bh2, bridges


def bicyclic_parent_name(core) -> str:
    """The 'bicyclo[x.y.z]alkane' parent name for `core`'s bridge lengths
    (P-23.2.2) — shared with `_von_baeyer_heteroatom.py`, since a skeletal
    heteroatom doesn't change the bridge-length count that determines this."""
    _, _, bridges = core
    lengths_desc = sorted((len(bridge) for bridge in bridges), reverse=True)
    total_atoms = sum(lengths_desc) + 2
    return f"bicyclo[{'.'.join(str(n) for n in lengths_desc)}]{alkane_name(total_atoms)}"


def iter_bicyclic_numberings(core):
    """Yield every von Baeyer-valid numbering (P-23.2.3) of `core` as a full
    atom-index order (position 0 -> locant 1, etc.) — every choice of
    starting bridgehead and, among bridges tied in length, their order,
    consistent with the fixed 'longer segment before shorter, main bridge
    last' shape. Shared with `_von_baeyer_heteroatom.py` so a heteroatom's
    locant can be minimized over the same candidate numberings a plain
    hydrocarbon's substituents are."""
    bh1, bh2, bridges = core
    for start, other, oriented_bridges in (
        (bh1, bh2, bridges),
        (bh2, bh1, [list(reversed(bridge)) for bridge in bridges]),
    ):
        for perm in permutations(range(3)):
            ordered = [oriented_bridges[i] for i in perm]
            if not (len(ordered[0]) >= len(ordered[1]) >= len(ordered[2])):
                continue
            main_ring_first, main_ring_second, main_bridge = ordered
            yield [start] + main_ring_first + [other] + list(reversed(main_ring_second)) + main_bridge


def _candidate_key(parent, substituents, heteroatom_locant=None, nondetachable_prefix=""):
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


def name_bicycloalkane(mol, core) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated bicyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    parent = bicyclic_parent_name(core)

    best_key = None
    best_name = None
    for full_order in iter_bicyclic_numberings(core):
        substituents = _substituents_for_ring(graph, full_order, halogens)
        key = _candidate_key(parent, substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, key[-1]

    return best_name
