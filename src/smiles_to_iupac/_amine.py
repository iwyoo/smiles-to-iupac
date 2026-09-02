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

This module accepts a *primary* amine (-NH2) attached to a non-aromatic
carbon, and (P-66.4/general N-substituent-prefix nomenclature, mirroring
`_amide.py`'s existing N-substituent handling) an acyclic secondary/tertiary
amine whose extra N-substituent(s) are simple unbranched, saturated alkyl
chains: the N-linked carbon starting the largest carbon skeleton becomes the
parent chain (suffixed '-amine' as usual), and each other N-linked chain is
cited as an 'N-'/'N,N-' substituent prefix, e.g. 'N-ethylethanamine'
(diethylamine, PubChem-verified) and 'N,N-dimethylmethanamine'
(trimethylamine, PubChem-verified). It otherwise mirrors `_alcohol.py`'s
scope restrictions:

Explicitly out of scope (raise `UnsupportedStructure`):
- A secondary/tertiary amine nitrogen with a branched, unsaturated, or
  ring-bearing N-substituent, or more than two N-substituents (mirrors
  `_amide.py`'s own N-substituent restrictions exactly).
- A secondary/tertiary amine nitrogen on or attached to a ring, or
  coexisting with another amine nitrogen elsewhere in the molecule (a
  diamine where one nitrogen is secondary/tertiary) — both deferred as
  separate, larger extensions.
- Any other heteroatom (O, S, ...) — including molecules that would also
  need a senior characteristic group (Table 3.3); this module rejects those
  outright rather than attempting suffix-vs-suffix seniority.
- -NH2/-NHR/-NR2 on an aromatic ring (aniline-type) — a separate module's
  territory.
- -NH2 on a von Baeyer polycyclic or spiro skeleton — deferred, same as
  `_alcohol.py`.
- An amine nitrogen on a carbon that is also part of a C=C/C#C bond (an
  enamine) — scoped out for the same reason `_alcohol.py` scopes out enols.

P-91.3/P-92: a molecule with one or
more *specified* tetrahedral stereocenters -- every one on the principal
chain/ring itself, no unspecified one alongside them, and no C=C/C#N
double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix, ascending
locant order, same pattern as `_thiol.py`/`_sulfonic_acid.py` (chain and
ring both). For a secondary/tertiary amine the stereodescriptor sits
outermost, ahead of the N-alkyl prefix, same as `_amide.py` (a
stereodescriptor always sits at the very front of the complete name,
P-91.3). The amine nitrogen itself is never a potential stereocenter
(pyramidal inversion), confirmed via RDKit `FindPotentialStereo` on
primary/secondary/tertiary examples alike.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    halogen_substituents,
    linear_branch,
    lowest_locant_set,
    non_single_bonds,
    path_between,
    specified_stereocenters,
)
from ._numerals import alkane_name, alkyl_name, multiplying_prefix, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _validate_and_collect_amines(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (amines, n_carbons_by_nitrogen): the set of amine-nitrogen
    atom indices, and a dict mapping each to a tuple of its 1-3 carbon
    neighbor indices (the -NH2 carbon for a primary amine; the parent-chain
    carbon plus 0-2 N-substituent carbons for a secondary/tertiary one)."""
    amines = set()
    n_carbons_by_nitrogen = {}
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an amine nitrogen (P-33.1) and "
                "halogen substituents (P-35.2.1) are not supported yet"
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
            neighbors = list(atom.GetNeighbors())
            if not neighbors or len(neighbors) > 3:
                raise UnsupportedStructure(
                    "a nitrogen with zero or more than three substituents "
                    "is not a valid amine nitrogen"
                )
            if any(n.GetAtomicNum() != 6 for n in neighbors):
                raise UnsupportedStructure(
                    "an amine nitrogen bonded to anything other than "
                    "carbon (e.g. another nitrogen) is out of scope for "
                    "this module"
                )
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() != 1.0:
                    raise UnsupportedStructure(
                        "a nitrogen double- or triple-bonded to carbon "
                        "(imine, nitrile) is not a plain amine, which this "
                        "module does not attempt to disambiguate"
                    )
            amines.add(atom.GetIdx())
            n_carbons_by_nitrogen[atom.GetIdx()] = tuple(n.GetIdx() for n in neighbors)
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
        raise UnsupportedStructure("no amine group found; this module only handles amines")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return amines, n_carbons_by_nitrogen


def _reject_enamine_carbons(graph, amines, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for n_idx in amines:
        for carbon in graph[n_idx]:
            if carbon in unsaturated_atoms:
                raise UnsupportedStructure(
                    "an amine nitrogen on a carbon that is also part of a "
                    "C=C/C#C bond (an enamine-type structure) is out of "
                    "scope for this module"
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
    """Each nitrogen's chain-side locant. `graph` may be the full-molecule
    adjacency (a secondary/tertiary nitrogen then has other, non-chain
    carbon neighbors too -- e.g. its N-substituents, deliberately excluded
    from the chain graph being tested), so a nitrogen only counts if exactly
    one of its carbon neighbors is part of the current candidate chain."""
    locants = []
    for n in amines:
        carbons_in_chain = [c for c in graph[n] if c in position_of]
        if len(carbons_in_chain) != 1:
            return None
        locants.append(position_of[carbons_in_chain[0]])
    return locants


def _substituents_for_chain(graph, chain, halogens, amines):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in amines]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _best_chain_name(carbon_graph, graph, halogens, amines, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch is out of scope). Returns (name, position_of) so a
    caller wrapping an N-alkyl prefix around `name` can still place a
    stereo prefix outermost using the winning chain's own locants."""
    chains = _longest_chains(carbon_graph)
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _amine_locants(position_of, amines, graph) is None:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        if stereo is not None and any(atom not in position_of for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every amine-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more -NH2 groups (P-44.1.1), an -NH2 expressed as "
            "an 'amino' substituent prefix, or a stereocenter on a "
            "substituent branch (P-92), is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            amine_locants = _amine_locants(position_of, amines, graph)
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, amines)
            key, name = _candidate_key(chain_length, amine_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of
    return best_name, best_position_of


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def _format_n_prefix(n_names):
    """Assemble the 'N-'/'N,N-'/'N,N,N-' prefix for a list of N-substituent
    names (one per "other" N-linked chain, duplicates included), grouping
    identical names under one multiplied prefix -- e.g. ['methyl',
    'methyl'] -> 'N,N-dimethyl' (trimethylamine), ['ethyl', 'methyl'] ->
    'N-ethyl-N-methyl' (asymmetric), ['methyl', 'methyl', 'methyl'] ->
    'N,N,N-trimethyl' (the quaternary tetramethylammonium cation's PIN
    'N,N,N-trimethylmethanaminium', P-73.1.1.1 -- confirmed against the
    primary source's own worked example, which explicitly rejects the
    'tetramethylazanium' alternative as non-PIN)."""
    counts = {}
    for name in n_names:
        counts[name] = counts.get(name, 0) + 1
    parts = []
    for name in sorted(counts):
        count = counts[name]
        if count == 1:
            parts.append(f"N-{name}")
        else:
            locants = ",".join(["N"] * count)
            parts.append(f"{locants}-{multiplying_prefix(count)}{name}")
    return "-".join(parts)


def _name_acyclic_secondary_tertiary_amine(mol, n_idx, n_carbons, bonds, stereo=None):
    """Name a secondary/tertiary amine: the N-linked carbon starting the
    largest carbon skeleton becomes the parent chain (suffixed '-amine' via
    `_best_chain_name`, same as a primary amine), and each other N-linked
    chain -- which must be a simple unbranched, saturated alkyl (mirrors
    `_amide.py`'s own N-substituent restriction) -- is cited as an
    'N-'/'N,N-' prefix (P-66.4), e.g. 'N-ethylethanamine' (diethylamine).

    Also reused directly by `_ammonium.py` for a quaternary ammonium
    cation's 3 "other" N-substituents (no charge/degree assumption is made
    here beyond the `n_carbons` tuple passed in -- `graph`/`full_carbon_graph`
    lookups are agnostic to the nitrogen's own formal charge)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    full_carbon_graph = carbon_adjacency(mol)

    # N itself isn't in the carbon-only graph, so the N-linked carbons never
    # touch each other there -- each already starts its own disjoint
    # component, letting the largest one be isolated as the parent chain's
    # own graph before any chain search runs.
    components = {c: _component_subgraph(full_carbon_graph, c) for c in n_carbons}
    parent_root = max(n_carbons, key=lambda c: len(components[c]))
    other_roots = [c for c in n_carbons if c != parent_root]

    n_names = []
    for other in other_roots:
        other_length = linear_branch(full_carbon_graph, other, None)
        if other_length is None:
            raise UnsupportedStructure("a branched N-substituent is not supported yet")
        other_atoms = set(components[other])
        if any(b[0] in other_atoms or b[1] in other_atoms for b in non_single_bonds(mol)):
            raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
        if any(neighbor in halogens for atom in other_atoms for neighbor in graph[atom]):
            raise UnsupportedStructure(
                "a halogen-substituted N-substituent is not supported yet "
                "-- it would be silently dropped, since only its carbon "
                "chain length is currently used to name it"
            )
        n_names.append(alkyl_name(other_length))

    excluded_atoms = set()
    for other in other_roots:
        excluded_atoms |= set(components[other])
    parent_carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in excluded_atoms}

    if any(neighbor in halogens for atom in parent_carbon_graph for neighbor in graph[atom]):
        raise UnsupportedStructure(
            "a secondary/tertiary amine coexisting with a halogen "
            "substituent on the parent chain is not supported yet -- "
            "P-14.5.2's alphanumerical interleaving of the 'N-' prefix "
            "with other substituent prefixes (e.g. PubChem's "
            "'2-chloro-N-ethylethanamine') is not implemented; this module "
            "would otherwise always cite the 'N-' prefix first"
        )

    best_name, best_position_of = _best_chain_name(parent_carbon_graph, graph, halogens, {n_idx}, bonds, stereo)

    n_prefix = _format_n_prefix(n_names)
    separator = "-" if best_name[0].isdigit() else ""
    best_name = f"{n_prefix}{separator}{best_name}"

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_acyclic_amine(mol, amines, n_carbons_by_nitrogen, bonds, stereo=None):
    if len(amines) == 1:
        (n_idx,) = amines
        n_carbons = n_carbons_by_nitrogen[n_idx]
        if len(n_carbons) > 1:
            return _name_acyclic_secondary_tertiary_amine(mol, n_idx, n_carbons, bonds, stereo)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    best_name, best_position_of = _best_chain_name(carbon_adjacency(mol), graph, halogens, amines, bonds, stereo)
    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
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


def _name_cyclic_amine(mol, amines, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must lie on
    the ring itself (P-92: a stereocenter on a substituent branch is out
    of scope, mirroring `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid`),
    and the winning ring numbering's own locants for those atoms are used
    to format a "(<locant><R/S>,...)-" prefix onto the name."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = _ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the ring "
            "itself is not supported yet (see P-92)"
        )

    best_key = None
    best_name = None
    best_position_of = None
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
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def name_amine(mol) -> str:
    amines, n_carbons_by_nitrogen = _validate_and_collect_amines(mol)
    stereo = specified_stereocenters(mol)
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
    has_secondary_or_tertiary = any(len(n_carbons_by_nitrogen[n]) > 1 for n in amines)
    if has_secondary_or_tertiary and num_rings > 0:
        raise UnsupportedStructure(
            "a secondary/tertiary amine nitrogen on or attached to a ring "
            "is out of scope for this module"
        )
    if has_secondary_or_tertiary and len(amines) > 1:
        raise UnsupportedStructure(
            "more than one amine nitrogen where at least one is "
            "secondary/tertiary is out of scope for this module"
        )
    if num_rings == 0:
        return _name_acyclic_amine(mol, amines, n_carbons_by_nitrogen, bonds, stereo)
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        return _name_cyclic_amine(mol, amines, stereo)
    raise UnsupportedStructure(
        "polycyclic and spiro amines are not supported yet (P-23/P-24/P-25 "
        "numbering integration with a suffix group is future work)"
    )
