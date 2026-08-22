"""Naming of simple monocyclic saturated hydrocarbons (cycloalkanes) whose ring
substituents are themselves unbranched, per the IUPAC 2013 Recommendations
("the Blue Book"):

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

Fused, bridged, and spiro ring systems, as well as branched ("compound")
substituents, are out of scope for this module and raise NotImplementedError.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    linear_branch,
    lowest_locant_set,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, alkyl_name, numerical_term


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


def _substituents_for_ring(graph, ring_order):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set]
        if not branch_roots:
            continue
        lengths = []
        for root in branch_roots:
            length = linear_branch(graph, root, atom)
            if length is None:
                return None
            lengths.append(length)
        substituents[position] = lengths
    return substituents


def _name_from_substituents(ring_size, substituents_by_name):
    parent = "cyclo" + alkane_name(ring_size)
    total_count = sum(len(v) for v in substituents_by_name.values())
    if total_count == 0:
        return parent
    if total_count == 1:
        # P-14.3.3: the locant is not essential on an otherwise unsubstituted ring.
        (name,) = substituents_by_name
        return f"{name}{parent}"
    prefixes = []
    for name in sorted(substituents_by_name):
        locants = sorted(substituents_by_name[name])
        multiplier = numerical_term(len(locants)) if len(locants) > 1 else ""
        prefixes.append(f"{','.join(str(loc) for loc in locants)}-{multiplier}{name}")
    return "-".join(prefixes) + parent


def _candidate_key(ring_size, substituents):
    """Sort key implementing P-45.2.2/P-45.2.3, most-preferred first (the
    substituent count is fixed for a given ring, so unlike the acyclic case
    there is no P-45.2.1 dimension to break ties on)."""
    substituents_by_name = {}
    for position, lengths in substituents.items():
        for length in lengths:
            substituents_by_name.setdefault(alkyl_name(length), []).append(position)

    locant_set = lowest_locant_set(
        loc for locants in substituents_by_name.values() for loc in locants
    )
    citation_locants = tuple(
        loc
        for name in sorted(substituents_by_name)
        for loc in sorted(substituents_by_name[name])
    )
    name = _name_from_substituents(ring_size, substituents_by_name)
    return locant_set, citation_locants, name


def name_cycloalkane(mol) -> str:
    validate_atoms_and_bonds(mol)
    if any(bond.GetBondTypeAsDouble() != 1.0 for bond in mol.GetBonds()):
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
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = _ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            substituents = _substituents_for_ring(graph, candidate)
            if substituents is None:
                continue
            key = _candidate_key(ring_size, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    if best_name is None:
        raise UnsupportedStructure(
            "branched (compound) substituent groups are not supported yet (see P-29.4)"
        )
    return best_name
