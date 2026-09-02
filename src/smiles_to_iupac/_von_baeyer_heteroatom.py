"""Skeletal replacement ('a' prefix) naming of a saturated von Baeyer
bicyclic (`_bicyclic.py`) or polycyclic, ring_count>=3 (`_polycyclic.py`)
ring system containing one or more ring heteroatoms, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-23.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  skeletal replacement ('a') nomenclature names a heteroatom in an
  otherwise-carbon von Baeyer ring system by citing a nondetachable
  prefix ('oxa', 'thia', 'aza', ...) with the locant of the replaced
  skeletal atom, directly in front of the 'bicyclo[x.y.z]alkane' (or
  'tricyclo[...]'/'tetracyclo[...]'/'pentacyclo[...]') parent (e.g.
  '7-oxabicyclo[2.2.1]heptane', '1-azabicyclo[2.2.2]octane' —
  quinuclidine). The bridge-length brackets and the ring's numbering
  (P-23.2.3) are unaffected by which skeletal atom is a heteroatom; only
  which of the numbering choices tied on that bracket is preferred
  changes: the heteroatom(s) must get the lowest locant set available
  among them, ahead of substituent locants (mirrors the substituent-locant
  tie-break already used for a plain hydrocarbon ring). For ring_count>=3
  this is *not* every numbering `_polycyclic.iter_polycyclic_candidates`
  yields, unlike the bicyclic case -- see
  `name_von_baeyer_heteroatom_polycyclic` below for why the topology
  choice (which candidates share the same bracket string) must still be
  settled first, by `iter_polycyclic_candidates`'s own `outer_key`.
- Two or more skeletal heteroatoms of the *same* element: P-14.2.1's
  ordinary multiplying prefix ('di', 'tri', ...) attaches directly to the
  'a'-term ('dioxa', 'triaza', ...), and every heteroatom's own locant is
  cited (ascending, comma-separated) before it -- confirmed against
  PubChem's own computed IUPACName for a from-scratch-built two-oxygen
  bicyclo[2.2.1]heptane, '2,7-dioxabicyclo[2.2.1]heptane' (C5H8O2,
  cross-checked further via ConnectivitySMILES). `_bicyclic._candidate_key`/
  `_polycyclic._candidate_key`'s existing `heteroatom_locant` parameter
  already ranks candidates correctly whether it's handed a single int or,
  as here, a locant tuple (Python's own tuple ordering does the "lowest
  locant set" comparison for free) -- no signature change was needed.
- Exactly two skeletal heteroatoms of *different* elements, one each:
  Table 2.8's element seniority order (P-23.2.1) is
  O > S > N among the three elements this module supports -- confirmed
  against PubChem's own computed IUPACName for three from-scratch-built
  bicyclo[3.2.1]octane variants sharing one skeleton
  ('8-oxa-3-azabicyclo[3.2.1]octane', '8-oxa-3-thiabicyclo[3.2.1]octane',
  '3-thia-8-azabicyclo[3.2.1]octane') plus a fourth, symmetric
  bicyclo[2.2.1]heptane built with O and N in mirror-image bridge
  positions ('2-oxa-6-azabicyclo[2.2.1]heptane') that isolates the tie-
  break: when the overall heteroatom locant *set* is tied between two
  numbering choices, the more senior element (not the lower atom index)
  gets the lower locant. The nondetachable prefix cites each element's
  own term ('oxa', 'thia', 'aza') with its locant, in seniority order
  (most senior first) regardless of which locant is numerically lower --
  e.g. '8-oxa-3-aza...', oxa first even though 3 < 8.

Explicitly out of scope (raise `UnsupportedStructure`):
- Three or more skeletal heteroatoms when not all the same element, or
  two heteroatoms of different elements where either element repeats
  (e.g. one O + two N) -- Table 2.8 seniority ordering above only covers
  exactly one heteroatom of each of two different elements.
- A heteroatom outside the ring core (e.g. a substituent -OH/-NH2
  alongside a plain hydrocarbon bicyclic — a different module's
  territory).
- Any heteroatom other than O, N, or S.
- Multiple heteroatoms in a ring_count>=3 polycyclic system (the
  ring_count>=3 case below still only accepts exactly one) -- this
  module's multi-heteroatom axis is bicyclic-only, and the mixed-element
  extension inherits that same bicyclic-only limit.
- Hexacyclic (ring_count=6) or larger polycyclic rings -- `_polycyclic.py`
  itself doesn't support these yet, independent of the heteroatom question.
- Unsaturation, charged/isotopic atoms, or anything else
  `_bicyclic.name_bicycloalkane`/`_polycyclic.name_polycycloalkane`'s own
  validation already rejects for the all-carbon case.
"""

from ._bicyclic import _candidate_key, bicyclic_parent_name, iter_bicyclic_numberings
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    multiplied_word,
    non_single_bonds,
)
from ._cyclic import _substituents_for_ring
from ._polycyclic import _candidate_key as _polycyclic_candidate_key, iter_polycyclic_candidates

_HETEROATOM_PREFIXES = {8: "oxa", 7: "aza", 16: "thia"}
_ALLOWED_ATOMIC_NUMS = {6, *_HETEROATOM_PREFIXES, *HALOGEN_PREFIXES}
# Table 2.8 (P-23.2.1) element seniority, restricted to the three elements
# this module supports: O > S > N (lower value = more senior).
_ELEMENT_SENIORITY = {8: 0, 16: 1, 7: 2}


def _validate_atoms(mol):
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a single skeletal O/N/S (P-23.2.1) "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num in HALOGEN_PREFIXES and atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1); "
                "polyvalent or bridging halogen structures are not supported"
            )


def _ring_heteroatoms(mol, core):
    bh1, bh2, bridges = core
    core_atoms = {bh1, bh2} | {atom for bridge in bridges for atom in bridge}
    return [atom for atom in core_atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() in _HETEROATOM_PREFIXES]


def has_single_ring_heteroatom_shape(mol, core) -> bool:
    """True iff `core` (a `_bicyclic.find_bicyclic_core` result) has exactly
    one O/N/S skeletal ring atom, the shape this module accepts."""
    return len(_ring_heteroatoms(mol, core)) == 1


def name_von_baeyer_heteroatom(mol, core) -> str:
    _validate_atoms(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated bicyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    ring_heteroatoms = _ring_heteroatoms(mol, core)
    all_heteroatoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    if len(ring_heteroatoms) != 1 or set(all_heteroatoms) != set(ring_heteroatoms):
        raise UnsupportedStructure(
            "exactly one skeletal ring heteroatom is supported (P-23.2.1's "
            "'a'-prefix ordering for two or more heteroatoms, and any "
            "heteroatom outside the ring skeleton, are out of scope here)"
        )
    (heteroatom_idx,) = ring_heteroatoms
    a_prefix = _HETEROATOM_PREFIXES[mol.GetAtomWithIdx(heteroatom_idx).GetAtomicNum()]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    parent = bicyclic_parent_name(core)

    best_key = None
    for full_order in iter_bicyclic_numberings(core):
        heteroatom_locant = full_order.index(heteroatom_idx) + 1
        substituents = _substituents_for_ring(graph, full_order, halogens)
        key = _candidate_key(
            parent,
            substituents,
            heteroatom_locant=heteroatom_locant,
            nondetachable_prefix=f"{heteroatom_locant}-{a_prefix}",
        )
        if best_key is None or key < best_key:
            best_key = key

    return best_key[-1]


def has_multi_ring_heteroatom_shape(mol, core) -> bool:
    """True iff `core` (a `_bicyclic.find_bicyclic_core` result) has two or
    more skeletal ring heteroatoms, all the same element -- the shape
    `name_von_baeyer_heteroatom_multi` accepts. A mixed-element ring (out of
    scope, Table 2.8 seniority not implemented) is deliberately excluded
    here rather than left to `name_von_baeyer_heteroatom_multi` to reject,
    so it falls through to whatever other module (if any) can name it."""
    ring_heteroatoms = _ring_heteroatoms(mol, core)
    if len(ring_heteroatoms) < 2:
        return False
    elements = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in ring_heteroatoms}
    return len(elements) == 1


def name_von_baeyer_heteroatom_multi(mol, core) -> str:
    _validate_atoms(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated bicyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    ring_heteroatoms = _ring_heteroatoms(mol, core)
    all_heteroatoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    elements = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in ring_heteroatoms}
    if len(ring_heteroatoms) < 2 or len(elements) != 1 or set(all_heteroatoms) != set(ring_heteroatoms):
        raise UnsupportedStructure(
            "two or more skeletal ring heteroatoms are only supported when "
            "all of the same element (P-23.2.1's 'a'-prefix ordering for "
            "mixed heteroatom kinds is out of scope here), and any "
            "heteroatom outside the ring skeleton is also out of scope"
        )
    (element,) = elements
    a_prefix = _HETEROATOM_PREFIXES[element]
    multiplied_a_prefix = multiplied_word(len(ring_heteroatoms), a_prefix)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    parent = bicyclic_parent_name(core)

    best_key = None
    for full_order in iter_bicyclic_numberings(core):
        heteroatom_locants = tuple(sorted(full_order.index(h) + 1 for h in ring_heteroatoms))
        substituents = _substituents_for_ring(graph, full_order, halogens)
        locant_citation = ",".join(str(loc) for loc in heteroatom_locants)
        key = _candidate_key(
            parent,
            substituents,
            heteroatom_locant=heteroatom_locants,
            nondetachable_prefix=f"{locant_citation}-{multiplied_a_prefix}",
        )
        if best_key is None or key < best_key:
            best_key = key

    return best_key[-1]


def has_mixed_element_heteroatom_shape(mol, core) -> bool:
    """True iff `core` (a `_bicyclic.find_bicyclic_core` result) has exactly
    two skeletal ring heteroatoms of two different elements, one each -- the
    shape `name_von_baeyer_heteroatom_mixed` accepts. Three or more
    heteroatoms, or two different elements where one repeats, are
    deliberately excluded (Table 2.8 seniority above only resolves the
    one-each case)."""
    ring_heteroatoms = _ring_heteroatoms(mol, core)
    if len(ring_heteroatoms) != 2:
        return False
    elements = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in ring_heteroatoms}
    return len(elements) == 2


def name_von_baeyer_heteroatom_mixed(mol, core) -> str:
    _validate_atoms(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated bicyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    ring_heteroatoms = _ring_heteroatoms(mol, core)
    all_heteroatoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    elements = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in ring_heteroatoms}
    if len(ring_heteroatoms) != 2 or len(elements) != 2 or set(all_heteroatoms) != set(ring_heteroatoms):
        raise UnsupportedStructure(
            "exactly one skeletal ring heteroatom of each of two different "
            "elements is supported (Table 2.8 seniority ordering for three "
            "or more heteroatoms, or a repeated element mixed with another, "
            "is out of scope here), and any heteroatom outside the ring "
            "skeleton is also out of scope"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    parent = bicyclic_parent_name(core)

    best_key = None
    for full_order in iter_bicyclic_numberings(core):
        by_seniority = sorted(
            ((mol.GetAtomWithIdx(a).GetAtomicNum(), full_order.index(a) + 1) for a in ring_heteroatoms),
            key=lambda pair: _ELEMENT_SENIORITY[pair[0]],
        )
        locant_set = tuple(sorted(loc for _, loc in by_seniority))
        heteroatom_key = (locant_set, tuple(loc for _, loc in by_seniority))
        nondetachable_prefix = "-".join(f"{loc}-{_HETEROATOM_PREFIXES[elem]}" for elem, loc in by_seniority)
        substituents = _substituents_for_ring(graph, full_order, halogens)
        key = _candidate_key(
            parent,
            substituents,
            heteroatom_locant=heteroatom_key,
            nondetachable_prefix=nondetachable_prefix,
        )
        if best_key is None or key < best_key:
            best_key = key

    return best_key[-1]


def _ring_heteroatoms_polycyclic(mol, core):
    branch_atoms, bridges = core
    core_atoms = set(branch_atoms) | {atom for _, _, path in bridges for atom in path}
    return [atom for atom in core_atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() in _HETEROATOM_PREFIXES]


def has_single_ring_heteroatom_shape_polycyclic(mol, core) -> bool:
    """True iff `core` (a `_polycyclic.find_polycyclic_core` result) has
    exactly one O/N/S skeletal ring atom, the shape this module accepts."""
    return len(_ring_heteroatoms_polycyclic(mol, core)) == 1


def name_von_baeyer_heteroatom_polycyclic(mol, core, ring_count) -> str:
    _validate_atoms(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated polycyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    ring_heteroatoms = _ring_heteroatoms_polycyclic(mol, core)
    all_heteroatoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    if len(ring_heteroatoms) != 1 or set(all_heteroatoms) != set(ring_heteroatoms):
        raise UnsupportedStructure(
            "exactly one skeletal ring heteroatom is supported (P-23.2.1's "
            "'a'-prefix ordering for two or more heteroatoms, and any "
            "heteroatom outside the ring skeleton, are out of scope here)"
        )
    (heteroatom_idx,) = ring_heteroatoms
    a_prefix = _HETEROATOM_PREFIXES[mol.GetAtomWithIdx(heteroatom_idx).GetAtomicNum()]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    for full_order, parent, outer_key in iter_polycyclic_candidates(core, ring_count):
        heteroatom_locant = full_order.index(heteroatom_idx) + 1
        substituents = _substituents_for_ring(graph, full_order, halogens)
        key = outer_key + _polycyclic_candidate_key(
            parent,
            substituents,
            heteroatom_locant=heteroatom_locant,
            nondetachable_prefix=f"{heteroatom_locant}-{a_prefix}",
        )
        if best_key is None or key < best_key:
            best_key = key

    return best_key[-1]
