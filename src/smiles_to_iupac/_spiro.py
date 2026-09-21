"""Naming of monospiro saturated hydrocarbons (two carbocyclic rings sharing
exactly one atom), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-24.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  "Monospiro parent hydrides consisting of two saturated cycloalkane rings
  are named by placing the nondetachable prefix 'spiro' before the name of
  the unbranched acyclic hydrocarbon with the same total number of skeletal
  atoms. The number of skeletal atoms linked to the spiro atom in each ring
  is indicated by arabic numbers separated by a full stop, cited in
  ascending order and enclosed in square brackets ... Numbering starts in
  the smaller ring, if one is smaller, at a ring atom next to the spiro atom
  and proceeds first around that ring, then through the spiro atom and
  around the second ring." (e.g. spiro[4.5]decane, spiro[4.4]nonane).
- P-24.2.0: this rule covers only saturated spiro systems built from two
  monocyclic rings; unsaturated alicyclic spiro ring systems are P-31.1.5
  and fused/bridged (non-spiro) polycyclic systems are P-23, both out of
  scope here.
- P-14.4 / P-45.2 (Chapter P-1 / P-4): when the base numbering path leaves a
  choice (which of the spiro atom's two ring neighbors starts each ring;
  which same-size ring is numbered first), lowest locants go to
  substituents as a set, then in order of citation — the same tie-break
  machinery `_cyclic.py` uses for monocyclic rings.
- P-29.4 / P-46 (Chapter P-2, P-4): branched ("compound") substituent
  groups — see `_substituents.py`.
- P-35.2.1 (Chapter P-3): halogen substituents (fluoro, chloro, bromo, iodo)
  hang off a ring atom the same way any other substituent does; no
  carbon-only filtering is needed here, same as in `_cyclic.py`.

Fused, bridged, and polyspiro ring systems are out of scope for this module
and raise UnsupportedStructure.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    non_single_bonds,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, substituents_for_ring


def find_monospiro_atom(mol):
    """Return the spiro atom index if `mol` is exactly two rings sharing one
    atom and no bonds (P-24.2.1's scope), else None."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 2:
        return None
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()
    shared_atoms = set(atom_rings[0]) & set(atom_rings[1])
    shared_bonds = set(bond_rings[0]) & set(bond_rings[1])
    if len(shared_atoms) != 1 or shared_bonds:
        return None
    return next(iter(shared_atoms))


def _walk_ring_from_spiro(graph, ring_set, spiro, start):
    full_ring_set = ring_set | {spiro}
    order = [start]
    previous, current = spiro, start
    while True:
        next_atom = next(n for n in graph[current] if n in full_ring_set and n != previous)
        if next_atom == spiro:
            return order
        order.append(next_atom)
        previous, current = current, next_atom


def _candidate_key(parent, substituents, heteroatom_locant=None, suffix_locant=None, nondetachable_prefix=""):
    """`heteroatom_locant`: shared with `_spiro_heteroatom.py`'s skeletal
    replacement-heteroatom case. `suffix_locant`: the analogous rank for a
    characteristic-group suffix (e.g. a monospiro alcohol's -OH,
    `_alcohol.py`'s `_name_monospiro_alcohol`) instead -- mirrors
    `_bicyclic._candidate_key`'s identical pair of parameters; the two
    never coexist here."""
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


def iter_monospiro_numberings(mol, spiro_atom):
    """Yield (parent, full_order) for every P-24.2.1-valid monospiro
    numbering: the ring-size-fixed choice of which ring is numbered first
    (smaller ring, or either when the two rings are tied) crossed with each
    ring's direction of traversal from the spiro atom. Shared with
    `_spiro_heteroatom.py` so a heteroatom's locant can be minimized over
    the same candidate numberings a plain hydrocarbon's substituents are."""
    graph = adjacency(mol)
    atom_rings = mol.GetRingInfo().AtomRings()
    ring_x = [atom for atom in atom_rings[0] if atom != spiro_atom]
    ring_y = [atom for atom in atom_rings[1] if atom != spiro_atom]

    if len(ring_x) < len(ring_y):
        ring_pairs = [(ring_x, ring_y)]
    elif len(ring_x) > len(ring_y):
        ring_pairs = [(ring_y, ring_x)]
    else:
        ring_pairs = [(ring_x, ring_y), (ring_y, ring_x)]

    for ring1, ring2 in ring_pairs:
        a, b = len(ring1), len(ring2)
        parent = f"spiro[{a}.{b}]{alkane_name(a + b + 1)}"
        ring1_set, ring2_set = set(ring1), set(ring2)
        start1 = next(n for n in graph[spiro_atom] if n in ring1_set)
        start2 = next(n for n in graph[spiro_atom] if n in ring2_set)
        order1 = _walk_ring_from_spiro(graph, ring1_set, spiro_atom, start1)
        order2 = _walk_ring_from_spiro(graph, ring2_set, spiro_atom, start2)
        for dir1 in (order1, list(reversed(order1))):
            for dir2 in (order2, list(reversed(order2))):
                yield parent, dir1 + [spiro_atom] + dir2


def name_monospiro(mol, spiro_atom) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated spiro ring systems are not supported yet (see "
            "P-31.1.5, unsaturated alicyclic spiro ring systems)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    best_name = None
    for parent, full_order in iter_monospiro_numberings(mol, spiro_atom):
        substituents = substituents_for_ring(graph, full_order, halogens)
        key = _candidate_key(parent, substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, key[-1]

    return best_name
