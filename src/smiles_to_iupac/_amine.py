"""Naming of primary amines (the '-amine' suffix, -NH2) on acyclic saturated
or unsaturated carbon chains and on simple monocyclic saturated rings, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-33.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'amine'
  is the preselected suffix for -NH2, cited in a combined chain/suffix name
  the same way 'ol' is for -OH (`_alcohol.py`), e.g. 'methanamine',
  'ethanamine', 'cyclohexanamine'.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission preconditions
  used in `_alcohol.py` for -OH apply identically to -NH2 (mononuclear
  parent, or a homogeneous two-carbon chain with exactly one substituent).
- P-44.4.1.8 / P-45.2: suffix locants are minimized before 'ene'/'yne'
  locants, which are minimized before substituent-prefix locants — same
  ordering as `_alcohol.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  -NH2 suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.

This module only accepts a *primary* amine (-NH2, nitrogen of degree 1 with
exactly two H's) attached to a non-aromatic carbon. It otherwise mirrors
`_alcohol.py`'s scope restrictions exactly:

Explicitly out of scope (raise `UnsupportedStructure`):
- Secondary/tertiary amines (nitrogen bonded to more than one carbon), or any
  nitrogen not shaped like an isolated -NH2.
- Any other heteroatom (O, S, ...) — including molecules that would also
  need a senior characteristic group (Table 3.3); this module rejects those
  outright rather than attempting suffix-vs-suffix seniority.
- -NH2 on an aromatic ring (aniline-type) — a separate module's territory.
- -NH2 on a von Baeyer polycyclic or spiro skeleton — deferred, same as
  `_alcohol.py`.
- -NH2 on a carbon that is also part of a C=C/C#C bond (an enamine) — scoped
  out for the same reason `_alcohol.py` scopes out enols.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    path_between,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _validate_and_collect_amines(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return the set of amine-nitrogen atom indices."""
    amines = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a primary amine nitrogen (P-33.1) "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 7:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a nitrogen bonded to more than one heavy atom (a "
                    "secondary/tertiary amine) is out of scope; only an "
                    "isolated primary amine (-NH2) is supported (P-33.1)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "a nitrogen double- or triple-bonded to carbon (imine, "
                    "nitrile) is not a plain primary amine, which this "
                    "module does not attempt to disambiguate"
                )
            if atom.GetTotalNumHs() != 2:
                raise UnsupportedStructure(
                    "an -N< atom that isn't a simple primary amine (-NH2) "
                    "is out of scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a primary amine must be attached to a carbon atom")
            amines.add(atom.GetIdx())
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if not amines:
        raise UnsupportedStructure("no primary amine (-NH2) group found; this module only handles amines")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return amines


def _reject_enamine_carbons(graph, amines, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for n_idx in amines:
        (carbon,) = graph[n_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a primary amine on a carbon that is also part of a "
                "C=C/C#C bond (an enamine-type structure) is out of scope "
                "for this module"
            )


def _multiplied_word(count, base):
    """P-16.3.3: a multiplying prefix's terminal 'a' is elided before a
    suffix beginning with 'a' or 'o' (see `_common.py`'s `multiplied_word`
    docstring for the confirmed examples this mirrors)."""
    if count == 0:
        return ""
    if count == 1:
        return base
    prefix = numerical_term(count)
    if prefix.endswith("a") and base[:1] in "ao":
        prefix = prefix[:-1]
    return prefix + base


def _suffix_body(ene_locants, yne_locants, amine_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'amine' endings
    (e.g. '4-en-1-amine'), plus whether the stem's trailing 'e' should be
    elided at the stem/first-segment boundary (only relevant when there is
    no 'ene'/'yne', see `_name_from_substituents`)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(amine_locants), _multiplied_word(len(amine_locants), "amine")))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    elide_stem = words[0][0] in "aeiouy"
    return "-".join(parts), elide_stem


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, amine_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants (prefix or suffix)
        # are always '1' and never cited.
        amine_word = _multiplied_word(len(amine_locants), "amine")
        stem = alkane_name(1)
        if amine_word[0] in "aeiouy":
            stem = stem[:-1]
        return format_substituent_prefixes(grouped, omit_locants=True) + stem + amine_word

    if chain_length == 2 and not has_unsaturation and total_subs == 0 and len(amine_locants) == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly one
        # substituent (here, the sole -NH2) in total has only one possible
        # structure, so the locant is omittable, e.g. 'ethanamine'.
        return alkane_name(2)[:-1] + "amine"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, amine_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, amine_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first."""
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    amine_locant_set = lowest_locant_set(amine_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, amine_locants, ene_locants, yne_locants, grouped)
    return (
        (
            amine_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _longest_chains(graph):
    nodes = list(graph)
    distances = {}
    parents = {}
    for node in nodes:
        dist, parent = bfs(graph, node)
        distances[node] = dist
        parents[node] = parent

    diameter = max(d for dist in distances.values() for d in dist.values())
    chains = []
    seen = set()
    for u in nodes:
        for v, d in distances[u].items():
            if d == diameter and (v, u) not in seen:
                seen.add((u, v))
                chains.append(path_between(parents[u], u, v))
    return chains


def _bond_locant(chain, bond_atoms):
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def _bond_locants(chain, bonds):
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _bond_locant(chain, (a, b))
        if locant is None:
            return None
        (ene if order == _ENE_ORDER else yne).append(locant)
    return ene, yne


def _amine_locants(position_of, amines, graph):
    locants = []
    for n in amines:
        (carbon,) = graph[n]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, amines):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in amines]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_amine(mol, amines, bonds):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _amine_locants(position_of, amines, graph) is None:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every amine-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more -NH2 groups (P-44.1.1), or an -NH2 expressed as "
            "an 'amino' substituent prefix, is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            amine_locants = _amine_locants(position_of, amines, graph)
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, amines)
            key, name = _candidate_key(chain_length, amine_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def _ring_cycle(graph, ring_atoms):
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def _substituents_for_ring(graph, ring_order, halogens, amines):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in amines]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, amine_locants, grouped):
    parent = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    amine_word = _multiplied_word(len(amine_locants), "amine")
    elide = amine_word[0] in "aeiouy"
    stem = parent[:-1] if elide else parent

    if total_subs == 0 and len(amine_locants) == 1:
        # P-14.3.3: the sole substituent on an otherwise unsubstituted ring
        # has no locant to distinguish, e.g. 'cyclohexanamine'.
        return stem + amine_word

    prefix = format_substituent_prefixes(grouped)
    loc_str = ",".join(str(loc) for loc in sorted(amine_locants))
    return f"{prefix}{stem}-{loc_str}-{amine_word}"


def _ring_candidate_key(ring_size, amine_locants, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    amine_locant_set = lowest_locant_set(amine_locants)
    name = _ring_name_from_substituents(ring_size, amine_locants, grouped)
    return amine_locant_set, locant_set, citation_locants, name


def _name_cyclic_amine(mol, amines):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = _ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            amine_locants = _amine_locants(position_of, amines, graph)
            if amine_locants is None:
                raise UnsupportedStructure(
                    "a primary amine not on the ring itself (e.g. on a "
                    "substituent branch) is not supported yet"
                )
            substituents = _substituents_for_ring(graph, candidate, halogens, amines)
            key = _ring_candidate_key(ring_size, amine_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_amine(mol) -> str:
    amines = _validate_and_collect_amines(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enamine_carbons(graph, amines, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_amine(mol, amines, bonds)
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        return _name_cyclic_amine(mol, amines)
    raise UnsupportedStructure(
        "polycyclic and spiro amines are not supported yet (P-23/P-24/P-25 "
        "numbering integration with a suffix group is future work)"
    )
