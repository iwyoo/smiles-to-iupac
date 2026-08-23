"""Naming of ketones (the '-one' suffix, C=O with two carbon substituents) on
acyclic saturated or unsaturated carbon chains and on simple monocyclic
saturated rings, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-33.4, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'one' is the preselected suffix for a ketone carbonyl, cited in a combined
  chain/suffix name the same way 'ol'/'amine' are for -OH/-NH2
  (`_alcohol.py`/`_amine.py`), e.g. 'propan-2-one', 'cyclohexanone'. Table
  3.3 ranks 'one' senior to 'ol' and 'amine' but junior to the carboxylic
  acid/ester/amide/nitrile/aldehyde suffixes; a carbonyl carbon shaped like
  an aldehyde (only one carbon neighbor) or a carboxylic acid/ester/amide
  (a second oxygen on the same carbon) is rejected outright rather than
  silently named as if it were a ketone. A coexisting hydroxyl (-OH), being
  junior to 'one', is *not* rejected: it is cited as the 'hydroxy'
  substituent prefix instead (P-41), e.g. 'CC(=O)CCO' ->
  '4-hydroxybutan-2-one'.
- P-44.4.1.8 / P-45.2: suffix locants are minimized before 'ene'/'yne'
  locants, which are minimized before substituent-prefix locants — same
  ordering as `_alcohol.py`/`_amine.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  ketone suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.

Unlike -OH/-NH2, a ketone carbonyl carbon always has exactly two carbon
neighbors and no hydrogens, so it can never be the sole substituent on a
mononuclear (P-14.3.4.2(a)) or homogeneous two-carbon (P-14.3.4.2(b)) chain
— those locant-omission special cases from `_alcohol.py` don't apply here
and are simply absent below; even the seemingly unambiguous 'propan-2-one'
(acetone) still cites its locant, confirmed against PubChem.

This module otherwise mirrors `_alcohol.py`'s scope restrictions:

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that isn't a doubly-bonded, isolated carbonyl oxygen or a
  singly-bonded hydroxyl (ethers, and any other oxygen shape).
- A carbonyl carbon with fewer than two carbon neighbors (aldehyde) or an
  aromatic carbonyl carbon (aryl ketone) — a separate module's territory.
- Any other heteroatom (N, S, ...) — this module only resolves the 'one'
  vs. 'ol' seniority competition (Table 3.3) between a ketone carbonyl and a
  coexisting hydroxyl; only C, halogen, and ketone/hydroxyl oxygen atoms are
  accepted at all.
- A hydroxyl on a carbon that is also part of a C=C/C#C bond (an enol,
  tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
- -one on a von Baeyer polycyclic or spiro skeleton — deferred, same as
  `_alcohol.py`.

Unlike -OH/-NH2, a ketone carbon can never itself also be a C=C/C#C alkene
carbon (its two remaining bonds, after the C=O double bond, are already
committed to its two required carbon substituents — a ketone carbon with a
third, double-bonded C=C neighbor would have five bonds), so there is no
enol-analogous "enone" case to scope out here: an alpha,beta-unsaturated
ketone like 'CC(=O)C=CC' (the double bond adjacent to, but not on, the
carbonyl carbon) is a perfectly nameable 'pent-3-en-2-one' and is supported
via the same suffix-locant-priority mechanism as `_alcohol.py`'s
'pent-4-en-1-ol'.
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
_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _validate_and_collect_ketones(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (ketones, hydroxyls): the set of carbonyl-oxygen atom indices,
    and the set of any coexisting hydroxyl-oxygen atom indices. A hydroxyl is
    junior to 'one' in Table 3.3's suffix seniority order, so it is cited as
    the 'hydroxy' substituent prefix instead of competing for the suffix
    (P-41)."""
    ketones = set()
    hydroxyls = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a ketone carbonyl oxygen (P-33.4) "
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
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether) is out of scope; only an isolated ketone "
                    "carbonyl or hydroxyl is supported (Table 3.3, P-33.4)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a ketone/hydroxyl oxygen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() == 1.0:
                if atom.GetTotalNumHs() != 1:
                    raise UnsupportedStructure(
                        "an oxygen that isn't a carbonyl (=O) or hydroxyl "
                        "(-OH) is out of scope for this module"
                    )
                hydroxyls.add(atom.GetIdx())
                continue
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring (an aryl ketone) is out "
                    "of scope for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                raise UnsupportedStructure(
                    "a carbonyl carbon with fewer than two carbon neighbors "
                    "(an aldehyde or terminal carbonyl) is a more senior "
                    "characteristic group than a plain ketone (Table 3.3), "
                    "which this module does not attempt to disambiguate"
                )
            ketones.add(atom.GetIdx())
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
    if not ketones:
        raise UnsupportedStructure("no ketone (C=O) group found; this module only handles ketones")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return ketones, hydroxyls


def _multiplied_word(count, base):
    if count == 0:
        return ""
    if count == 1:
        return base
    return numerical_term(count) + base


def _suffix_body(ene_locants, yne_locants, one_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'one' endings
    (e.g. '4-en-1-one')."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(one_locants), _multiplied_word(len(one_locants), "one")))

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


def _name_from_substituents(chain_length, one_locants, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, one_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, one_locants, ene_locants, yne_locants, substituents):
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
    one_locant_set = lowest_locant_set(one_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, one_locants, ene_locants, yne_locants, grouped)
    return (
        (
            one_locant_set,
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


def _one_locants(position_of, ketones, graph):
    locants = []
    for o in ketones:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, ketones):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in ketones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_ketone(mol, ketones, hydroxyls, bonds):
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _one_locants(position_of, ketones, graph) is None:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every ketone-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more C=O groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ketones, graph)
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, ketones)
            key, name = _candidate_key(chain_length, one_locants, ene_locants, yne_locants, substituents)
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


def _substituents_for_ring(graph, ring_order, halogens, ketones):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in ketones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, one_locants, grouped):
    parent = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    one_word = _multiplied_word(len(one_locants), "one")
    elide = one_word[0] in "aeiouy"
    stem = parent[:-1] if elide else parent

    if total_subs == 0 and len(one_locants) == 1:
        # P-14.3.3: the sole substituent on an otherwise unsubstituted ring
        # has no locant to distinguish, e.g. 'cyclohexanone'.
        return stem + one_word

    prefix = format_substituent_prefixes(grouped)
    loc_str = ",".join(str(loc) for loc in sorted(one_locants))
    return f"{prefix}{stem}-{loc_str}-{one_word}"


def _ring_candidate_key(ring_size, one_locants, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    one_locant_set = lowest_locant_set(one_locants)
    name = _ring_name_from_substituents(ring_size, one_locants, grouped)
    return one_locant_set, locant_set, citation_locants, name


def _name_cyclic_ketone(mol, ketones, hydroxyls):
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
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
            one_locants = _one_locants(position_of, ketones, graph)
            if one_locants is None:
                raise UnsupportedStructure(
                    "a ketone not on the ring itself (e.g. on a substituent "
                    "branch) is not supported yet"
                )
            substituents = _substituents_for_ring(graph, candidate, halogens, ketones)
            key = _ring_candidate_key(ring_size, one_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_ketone(mol) -> str:
    ketones, hydroxyls = _validate_and_collect_ketones(mol)
    graph = adjacency(mol)
    # Exclude each C=O carbonyl bond itself: `non_single_bonds` reports it as
    # order 2.0 same as a C=C, but it isn't a chain 'ene' bond (one endpoint
    # is the ketone oxygen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in ketones and b[1] not in ketones]
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior ketone/aldehyde "
                "form and is out of scope for this module (P-31.1.4.2.4)"
            )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_ketone(mol, ketones, hydroxyls, bonds)
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        return _name_cyclic_ketone(mol, ketones, hydroxyls)
    raise UnsupportedStructure(
        "polycyclic and spiro ketones are not supported yet (P-23/P-24/P-25 "
        "numbering integration with a suffix group is future work)"
    )
