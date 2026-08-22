"""Naming of simple monocyclic saturated hydrocarbons (cycloalkanes), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-22.1.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): the name
  is formed by attaching the nondetachable prefix 'cyclo' to the name of the
  acyclic saturated hydrocarbon with the same number of carbon atoms.
- P-14.4 / P-45.2 (Chapter P-1 / P-4): the numbering direction and starting atom
  are chosen, as for a chain, to give the lowest locant set to substituents,
  then the lowest locants in their order of citation.
- P-14.3.3 (Chapter P-1): a locant is cited only when essential to the
  structure; for a single substituent on an otherwise unsubstituted ring, every
  ring atom is equivalent before substitution, so the locant is not essential
  and is omitted (e.g. 'methylcyclohexane', not '1-methylcyclohexane').
- P-29.4 / P-46 (Chapter P-2, P-4): branched ("compound") substituent groups,
  e.g. `(1-methylpropyl)` for a sec-butyl-like ring substituent — see
  `_substituents.py`.
- P-35.2.1 (Chapter P-3): halogen substituents (fluoro, chloro, bromo, iodo)
  hang off a ring atom the same way any other substituent does; a ring's own
  atom sequence needs no carbon-only filtering here since RDKit's ring
  perception (`GetRingInfo`) never includes a monovalent atom in a ring.

Fused, bridged, and spiro ring systems are out of scope for this module and
raise NotImplementedError.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch


def _ring_cycle(graph, ring_atoms):
    """Order a simple ring's atoms into a cyclic sequence by walking its bonds."""
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def _substituents_for_ring(graph, ring_order, halogens):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set]
        if not branch_roots:
            continue
        substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(ring_size, grouped):
    parent = "cyclo" + alkane_name(ring_size)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    if total_count == 0:
        return parent
    if total_count == 1:
        # P-14.3.3: the locant is not essential on an otherwise unsubstituted ring.
        (name,) = grouped
        display_name = f"({name})" if grouped[name]["compound"] else name
        return f"{display_name}{parent}"
    prefix = format_substituent_prefixes(grouped)
    return prefix + parent


def _candidate_key(ring_size, substituents):
    """Sort key implementing P-45.2.2/P-45.2.3, most-preferred first (the
    substituent count is fixed for a given ring, so unlike the acyclic case
    there is no P-45.2.1 dimension to break ties on)."""
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(ring_size, grouped)
    return locant_set, citation_locants, name


def name_cycloalkane(mol) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated rings are not supported yet (see P-31.1.3, "
            "cycloalkenes and cycloalkynes)"
        )
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        raise UnsupportedStructure(
            "only simple monocyclic ring systems are supported so far (see "
            "P-23/P-24/P-25 for polyalicyclic, spiro, and fused ring systems)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = _ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            substituents = _substituents_for_ring(graph, candidate, halogens)
            key = _candidate_key(ring_size, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    return best_name
