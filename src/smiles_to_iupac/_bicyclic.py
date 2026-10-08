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
- P-31.1.4.1/P-31.1.4.2 (`_name_bicyclic_unsaturated`): a bicyclic ring
  system bearing one or more C=C/C#C double/triple bonds still uses the
  identical 'bicyclo[x.y.z]' numbering above -- the 'ane' ending is simply
  replaced with 'ene'/'yne', and the bond locant(s) are folded into the
  existing numbering tie-break ahead of substituent locants, mirroring
  `_von_baeyer_heteroatom.py`'s identical heteroatom-locant priority.
  P-31.1.4.1's simple case (every bond's two atoms land on numerically
  adjacent locants under some valid numbering, no ring-wraparound the way
  a plain monocyclic ring allows) is cited as a plain locant, confirmed
  against the Blue Book's own worked example 'bicyclo[3.2.1]oct-2-ene
  (PIN)'. P-31.1.4.2's compound-locant case (no numbering makes a bond's
  atoms adjacent) cites the higher locant in parentheses instead (e.g.
  'bicyclo[2.2.1]hept-1(7)-ene', PubChem CID 53949421), ranked by its own
  3-step tie-break (`_common.von_baeyer_unsaturation_citations`).
  A true Kekule-aromatic-ring case (a literally-aromatic-flagged ring
  fused/bridged into the system) is a separate, further follow-up.

Tricyclic and higher polycyclic systems (P-23.2.5, P-23.2.6) and
heteroatom-containing bicyclics (P-23.3) are out of scope and raise
UnsupportedStructure.
"""

from itertools import permutations

from ._common import (
    ENE_BOND_ORDER,
    YNE_BOND_ORDER,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    kekulized_copy,
    lowest_locant_set,
    non_single_bonds,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
    specified_double_bond_stereo,
    von_baeyer_bond_stereo,
    von_baeyer_unsaturation_citations,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, substituents_for_ring
from ._unsaturated import _unsaturation_suffix_from_citations


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


def _candidate_key(parent, substituents, heteroatom_locant=None, suffix_locant=None, nondetachable_prefix=""):
    """`suffix_locant`: like `heteroatom_locant` (same tie-break rank,
    ahead of substituent locants) but for a characteristic-group suffix
    (e.g. a von Baeyer alcohol's -OH, `_alcohol.py`'s `_name_von_baeyer_
    alcohol`) instead of a skeletal replacement heteroatom -- the two
    never coexist in this codebase yet (a plain carbocyclic ring's own
    suffix vs. a skeletal-replacement heteroatom), so their relative
    priority when both are given is left unresolved."""
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    if prefix and nondetachable_prefix:
        prefix += "-"
    name = prefix + nondetachable_prefix + parent
    if heteroatom_locant is not None:
        return heteroatom_locant, locant_set, citation_locants, name
    if suffix_locant is not None:
        return suffix_locant, locant_set, citation_locants, name
    return locant_set, citation_locants, name


def _name_bicyclic_unsaturated(mol, core, bonds) -> str:
    """P-31.1.4.1's simple case (every double/triple bond's two atoms land
    on consecutive locants) and P-31.1.4.2's compound-locant case
    (otherwise, the higher locant cited in parentheses) both handled
    uniformly via `_common.von_baeyer_unsaturation_citations`: every
    candidate numbering (`iter_bicyclic_numberings`) always yields *some*
    valid citation now, ranked by P-31.1.4.2's own 3-step tie-break --
    (1) fewest compound locants, (2) lowest primary (parenthesized-number-
    ignoring) locant set, (3) lowest full locant set including
    parenthesized numbers -- folded into `_candidate_key`'s existing
    `suffix_locant` slot ahead of substituent locants."""
    invalid_orders = [order for _, _, order in bonds if order not in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if invalid_orders:
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    stem = bicyclic_parent_name(core)[:-3]
    bond_stereo = specified_double_bond_stereo(mol)

    best_key = None
    for full_order in iter_bicyclic_numberings(core):
        position = {atom: i + 1 for i, atom in enumerate(full_order)}
        z_locants, stereo_prefix = von_baeyer_bond_stereo(mol, position, bond_stereo)
        ene_citations, yne_citations, compound_count, primary_locants, full_locants = (
            von_baeyer_unsaturation_citations(position, bonds)
        )
        primary_locant_set = lowest_locant_set(primary_locants)
        full_locant_set = lowest_locant_set(full_locants)
        ene_locant_set = lowest_locant_set([locant for locant, _ in ene_citations])
        body, needs_stem_a = _unsaturation_suffix_from_citations(ene_citations, yne_citations)
        parent = stem + ("a" if needs_stem_a else "") + "-" + body
        substituents = substituents_for_ring(graph, full_order, halogens, mol=mol, unsaturated=True)
        key = _candidate_key(
            parent,
            substituents,
            suffix_locant=(compound_count, primary_locant_set, full_locant_set, ene_locant_set, z_locants),
        )
        key = key[:-1] + (stereo_prefix + key[-1],)
        if best_key is None or key < best_key:
            best_key = key

    return best_key[-1]


def name_bicycloalkane(mol, core) -> str:
    validate_atoms_and_bonds(mol)
    if any(a.GetIsAromatic() for a in mol.GetAtoms()) and bicyclic_parent_name(core).endswith(".0]" + bicyclic_parent_name(core).split("]")[-1]):
        raise UnsupportedStructure("an aromatic ortho-fused bicyclic ring system is named by fusion nomenclature, not von Baeyer")
    mol = kekulized_copy(mol)
    ring_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
    bonds = [b for b in non_single_bonds(mol) if b[0] in ring_atoms and b[1] in ring_atoms]
    if bonds:
        return _name_bicyclic_unsaturated(mol, core, bonds)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    parent = bicyclic_parent_name(core)

    best_key = None
    best_name = None
    for full_order in iter_bicyclic_numberings(core):
        substituents = substituents_for_ring(graph, full_order, halogens, mol=mol, unsaturated=True)
        key = _candidate_key(parent, substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, key[-1]

    return best_name
