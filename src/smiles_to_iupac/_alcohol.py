"""Naming of alcohols (the '-ol' suffix, -OH) on acyclic saturated or
unsaturated carbon chains and on simple monocyclic saturated rings, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-33.2.1, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'ol' is the preselected suffix for -OH, ranked 14th (out of 17) in Table
  3.3's seniority-for-citation-as-suffix order. Table 3.3 also lists several
  suffixes senior to 'ol' that this module must not be fooled into
  misreading as a plain alcohol: carboxylic acid/oic acid ('-CO-OH'), amide,
  nitrile, and in particular 'al' ('-CHO') and 'one' ('=O') — any oxygen
  double-bonded to a carbon (aldehyde, ketone, or the carbonyl half of a
  carboxylic acid) is rejected outright rather than silently named as if it
  were an alcohol.
- P-30 Introduction (Chapter P-3): "a characteristic group systematically
  introduced as a suffix attached to a parent hydride, for example ...
  ethanol" is the Blue Book's own example of this construction: the parent
  hydride name ('ethane') with its final 'e' elided before the vowel-initial
  suffix 'ol' -> 'ethanol'.
- P-44.1.1 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): the
  senior parent structure has the *maximum number* of principal
  characteristic groups (here, -OH) — this is applied before chain length
  (P-44.3.2) is even considered, e.g. 'pentane-1,4-diol' is preferred over a
  longer chain that would only capture one of two -OH groups on its own
  chain, the second -OH being pushed into a branch as a 'hydroxy' prefix.
  This module does NOT implement that full generality: it only accepts a
  candidate chain that is simultaneously (a) one of the longest carbon
  chains (P-44.3.2, as in `_acyclic.py`/`_unsaturated.py`) and (b) carries
  every -OH-bearing carbon in the molecule (mirroring `_unsaturated.py`'s
  existing, analogous restriction that every multiple bond must lie on a
  candidate longest chain). A molecule where the true PIN would require a
  shorter chain to capture more -OH groups, or would push a leftover -OH
  into a 'hydroxy' substituent prefix, is out of scope and raises
  `UnsupportedStructure`.
- P-44.4.1 (Chapter P-4): once the parent chain/ring is otherwise fixed, its
  *numbering* is chosen by criteria (a)-(l) in order; in particular
  criterion (h) (P-44.4.1.8, "the lower locant for an attached group
  expressed as a suffix") outranks criterion (j) (P-44.4.1.10, the 'ene'/
  'yne' locants) which in turn outranks substituent-prefix locants (P-45.2,
  applied after all of P-44.4.1 — see `_acyclic.py`/`_unsaturated.py`).
  Concretely: the -OH locant set is minimized first; only if that leaves a
  choice are 'ene'/'yne' locants minimized; only if that still leaves a
  choice are substituent-prefix locants minimized. Example confirmed against
  this exact ordering: 'pent-4-en-1-ol' (not 'pent-1-en-5-ol' /
  'pent-1-en-4-ol'-numbered-from-the-alkene-end) — the -OH gets locant 1
  even though this gives the double bond the *higher* available locant (4
  rather than 1).
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-alcohol name (e.g. 'pent-4-en-1-ol')
  reuses the same mechanics as `_unsaturated.py` (ending replaces 'ane'
  entirely, multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene'
  before 'yne', euphonic stem 'a' before a multiplied ending).
- P-14.3.4.2(a)/(b) (Chapter P-1, as already applied to halogens in
  `_acyclic.py`): a locant is omitted from a mononuclear (one-carbon) parent
  regardless of how many substituents/suffixes it carries ('methanol'), and
  from a homogeneous two-carbon chain bearing *exactly one* substituent in
  total, suffix or prefix ('ethanol'; 'ethane-1,2-diol' still needs locants
  because it has two -OH's, not one). This module applies the same
  precondition uniformly to prefix and suffix locants alike, since 'ethanol'
  is itself the Blue Book's own worked example of this omission applied to a
  suffix (see P-30 above). NOTE: for a two-carbon chain carrying exactly one
  halogen prefix *and* the -OH suffix together (two different substituents,
  so the omission precondition above is not met), whether the -OH locant
  must still be cited in the strict PIN is not fully resolved here — PubChem's
  auto-generated name for this exact case ('2-chloroethanol') omits it, but
  this module's tests deliberately avoid that specific ambiguous case and
  this module instead always cites both locants once the omission
  precondition fails (e.g. '3-chloropropan-1-ol', unambiguous since a
  three-carbon-or-longer chain never qualifies for omission at all), for
  consistency with the general locant-citation rule already used everywhere
  else in this project.
- P-35.2.1 (Chapter P-3): halogen substituents are prefix-only and coexist
  freely with the -OH suffix (they never compete for suffix status), reusing
  `halogen_substituents`/`format_substituent_prefixes` unchanged.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that is not an isolated, singly-bonded -OH with exactly one H
  (ethers, and any C=O — aldehyde, ketone, or the carbonyl of a carboxylic
  acid/amide/ester).
- Two or more *different* characteristic-group types (e.g. an alcohol and an
  amine) — not applicable here since only C, halogen, and -OH-shaped oxygen
  atoms are accepted at all; any other heteroatom (N, S, ...) is rejected.
- -OH on an aromatic ring (phenol-type) — a separate, in-progress module's
  territory.
- -OH on a von Baeyer polycyclic (bicyclic through pentacyclic) or spiro
  skeleton — deferred; those modules' internal numbering would need real
  integration work to prioritize a suffix locant correctly.
- -OH on a carbon that is also part of a C=C/C#C bond (an enol) — this is a
  further, deliberate narrowing beyond what full generality would allow (the
  Blue Book's numbering machinery *would* handle this case, combining suffix
  and 'ene'/'yne' locant priority the same way as any other case), but the
  keto-enol-tautomer-adjacent naming nuances were not independently verified
  here, so it is scoped out rather than guessed at.
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


def _validate_and_collect_hydroxyls(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return the set of hydroxyl-oxygen atom indices."""
    hydroxyls = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a hydroxyl oxygen (P-33.2.1) and "
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
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether) is out of scope; only an isolated hydroxyl "
                    "(-OH) is supported (Table 3.3, P-33.2.1)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "an oxygen double-bonded to carbon indicates a more "
                    "senior characteristic group (aldehyde, ketone, or a "
                    "carboxylic acid) than a plain alcohol, which this "
                    "module does not attempt to disambiguate (Table 3.3's "
                    "suffix seniority order, P-33.2.1)"
                )
            if atom.GetTotalNumHs() != 1:
                raise UnsupportedStructure(
                    "an -O- atom that isn't a simple hydroxyl (-OH) is out "
                    "of scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a hydroxyl must be attached to a carbon atom")
            hydroxyls.add(atom.GetIdx())
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
    if not hydroxyls:
        raise UnsupportedStructure("no hydroxyl (-OH) group found; this module only handles alcohols")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return hydroxyls


def _reject_enol_carbons(graph, hydroxyls, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o_idx in hydroxyls:
        (carbon,) = graph[o_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol-type structure) is out of scope for this module"
            )


def _multiplied_word(count, base):
    if count == 0:
        return ""
    if count == 1:
        return base
    return numerical_term(count) + base


def _suffix_body(ene_locants, yne_locants, oh_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'ol' endings
    (e.g. '4-en-1-ol'), plus whether the stem's trailing 'e' should be
    elided at the stem/first-segment boundary (only relevant when there is
    no 'ene'/'yne', see `_name_from_substituents`)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(oh_locants), _multiplied_word(len(oh_locants), "ol")))

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


def _name_from_substituents(chain_length, oh_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants (prefix or suffix)
        # are always '1' and never cited.
        ol_word = _multiplied_word(len(oh_locants), "ol")
        stem = alkane_name(1)
        if ol_word[0] in "aeiouy":
            stem = stem[:-1]
        return format_substituent_prefixes(grouped, omit_locants=True) + stem + ol_word

    if chain_length == 2 and not has_unsaturation and total_subs == 0 and len(oh_locants) == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly one
        # substituent (here, the sole -OH) in total has only one possible
        # structure, so the locant is omittable, e.g. 'ethanol (PIN)'.
        return alkane_name(2)[:-1] + "ol"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, oh_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, oh_locants, ene_locants, yne_locants, substituents):
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
    oh_locant_set = lowest_locant_set(oh_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, oh_locants, ene_locants, yne_locants, grouped)
    return (
        (
            oh_locant_set,
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


def _oh_locants(position_of, hydroxyls, graph):
    locants = []
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, hydroxyls):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in hydroxyls]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_alcohol(mol, hydroxyls, bonds):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _oh_locants(position_of, hydroxyls, graph) is None:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every hydroxyl-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more -OH groups (P-44.1.1), or an -OH expressed as a "
            "'hydroxy' substituent prefix, is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            oh_locants = _oh_locants(position_of, hydroxyls, graph)
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, hydroxyls)
            key, name = _candidate_key(chain_length, oh_locants, ene_locants, yne_locants, substituents)
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


def _substituents_for_ring(graph, ring_order, halogens, hydroxyls):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in hydroxyls]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, oh_locants, grouped):
    parent = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    ol_word = _multiplied_word(len(oh_locants), "ol")
    elide = ol_word[0] in "aeiouy"
    stem = parent[:-1] if elide else parent

    if total_subs == 0 and len(oh_locants) == 1:
        # P-14.3.3: the sole substituent on an otherwise unsubstituted ring
        # has no locant to distinguish, e.g. 'cyclohexanol'.
        return stem + ol_word

    prefix = format_substituent_prefixes(grouped)
    loc_str = ",".join(str(loc) for loc in sorted(oh_locants))
    return f"{prefix}{stem}-{loc_str}-{ol_word}"


def _ring_candidate_key(ring_size, oh_locants, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    oh_locant_set = lowest_locant_set(oh_locants)
    name = _ring_name_from_substituents(ring_size, oh_locants, grouped)
    return oh_locant_set, locant_set, citation_locants, name


def _name_cyclic_alcohol(mol, hydroxyls):
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
            oh_locants = _oh_locants(position_of, hydroxyls, graph)
            if oh_locants is None:
                raise UnsupportedStructure(
                    "a hydroxyl not on the ring itself (e.g. on a "
                    "substituent branch) is not supported yet"
                )
            substituents = _substituents_for_ring(graph, candidate, halogens, hydroxyls)
            key = _ring_candidate_key(ring_size, oh_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_alcohol(mol) -> str:
    hydroxyls = _validate_and_collect_hydroxyls(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enol_carbons(graph, hydroxyls, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_alcohol(mol, hydroxyls, bonds)
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        return _name_cyclic_alcohol(mol, hydroxyls)
    raise UnsupportedStructure(
        "polycyclic and spiro alcohols are not supported yet (P-23/P-24/"
        "P-25 numbering integration with a suffix group is future work)"
    )
