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


def name_monospiro(mol, spiro_atom) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated spiro ring systems are not supported yet (see "
            "P-31.1.5, unsaturated alicyclic spiro ring systems)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    atom_rings = mol.GetRingInfo().AtomRings()
    ring_x = [atom for atom in atom_rings[0] if atom != spiro_atom]
    ring_y = [atom for atom in atom_rings[1] if atom != spiro_atom]

    if len(ring_x) < len(ring_y):
        ring_pairs = [(ring_x, ring_y)]
    elif len(ring_x) > len(ring_y):
        ring_pairs = [(ring_y, ring_x)]
    else:
        ring_pairs = [(ring_x, ring_y), (ring_y, ring_x)]

    best_key = None
    best_name = None
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
                full_order = dir1 + [spiro_atom] + dir2
                substituents = _substituents_for_ring(graph, full_order, halogens)
                key = _candidate_key(parent, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]

    return best_name
