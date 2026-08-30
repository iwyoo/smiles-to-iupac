"""Naming of sulfonamides (the '-sulfonamide' suffix, -SO2NH2) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.3.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf), the
  nitrogen analogue of the sulfonic acid entry in the same table:
  'sulfonamide' is the preselected suffix for -S(=O)(=O)-NH2, mirroring
  `_sulfonic_acid.py`'s own '-SO3H' shape with the hydroxyl oxygen replaced
  by a primary amide nitrogen. This module doesn't implement acid/amide-vs-
  other-suffix seniority yet (see scope note below), mirroring
  `_sulfonic_acid.py`'s own deferral.
- Like 'sulfonic acid', 'sulfonamide' is cited as the parent hydride name
  followed directly by the suffix with no elision -- 'methane' +
  'sulfonamide' -> 'methanesulfonamide' (PubChem structure match).
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_sulfonic_acid.py` apply (mononuclear parent, or a homogeneous two-carbon
  chain with exactly one substituent in total), e.g. 'ethanesulfonamide'
  (PubChem structure match).
- P-44.4.1.8 / P-45.2: the -SO2NH2 locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as `_sulfonic_acid.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -SO2NH2 suffix. Confirmed via PubChem: 'propane-1-sulfonamide'
  (CCCS(=O)(=O)N), 'cyclohexanesulfonamide' (O=S(=O)(N)C1CCCCC1),
  '2-chloroethanesulfonamide' (ClCCS(=O)(=O)N).
- P-91.3/P-92 (`tasks/sulfonamide-stereocenter-naming.md`): a molecule with
  one or more *specified* tetrahedral stereocenters -- every one on the
  principal chain/ring itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
  ascending locant order, same pattern as `_sulfonic_acid.py` (chain and
  ring both). The sulfonamide's own sulfur is never itself a potential
  stereocenter (its two =O substituents are identical), confirmed via
  RDKit `FindPotentialStereo` on `CC(C)S(=O)(=O)N`.
- The sulfonamide nitrogen may carry zero, one, or two plain, unbranched,
  unsubstituted, saturated alkyl substituents, each cited as its own
  'N-'-prefixed substituent directly ahead of the parent stem, in
  alphabetical order, with a 'di' multiplying prefix (and a single shared
  'N,N-' pair) when both are identical -- exactly `_amide.py`'s own N-/
  N,N-disubstitution pattern, reusing `alkyl_name` directly (not
  `name_branch`, see the P-29 blocker) since both substituents are
  restricted to an unbranched chain. Confirmed via PubChem:
  'N-methylmethanesulfonamide' (CS(=O)(=O)NC), 'N,N-
  dimethylmethanesulfonamide' (CS(=O)(=O)N(C)C).

Scope, deliberately narrow, mirroring `_sulfonic_acid.py`'s own first pass
exactly: a single -SO2NH2 on an acyclic chain or on a single saturated
carbon ring (monocyclic), with no other heteroatom anywhere in the molecule
except the sulfonamide group's own oxygens/nitrogen (and any N-alkyl
substituent's carbons) -- acid/amide-vs-other Table 3.3 seniority
coexistence is future work. Explicitly out of scope (raise
`UnsupportedStructure`): a branched, unsaturated, or ring-bearing
N-substituent, polycyclic/spiro/unsaturated rings, a -SO2NH2 on a
substituent branch off an otherwise-unsubstituted ring, two or more
-SO2NH2 groups, and a sulfonamide on a carbon that is also part of a
C=C/C#C bond.
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
    ring_cycle,
    specified_stereocenters,
)
from ._numerals import alkane_name, alkyl_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _sulfonamide_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfonamide group: bonded to exactly one
    carbon, two double-bonded (terminal) oxygens, and one single-bonded
    amide nitrogen (0-2 carbon substituents besides the sulfur, otherwise
    terminal with the rest hydrogens)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(oxygens) != 2 or len(nitrogens) != 1:
            continue
        double_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0]
        if len(double_os) != 2 or any(o.GetDegree() != 1 for o in double_os):
            continue
        (nitrogen,) = nitrogens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        n_substituents = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        if len(n_substituents) > 2 or any(n.GetAtomicNum() != 6 for n in n_substituents):
            continue
        if nitrogen.GetTotalNumHs() != 2 - len(n_substituents):
            continue
        if any(bond.GetBondTypeAsDouble() != 1.0 for bond in nitrogen.GetBonds()):
            continue
        matches.append(atom)
    return matches


def has_sulfonamide_shape(mol) -> bool:
    return bool(_sulfonamide_sulfur_atoms(mol))


def _validate_and_collect_sulfonamides(mol):
    sulfur_atoms = _sulfonamide_sulfur_atoms(mol)
    if not sulfur_atoms:
        raise UnsupportedStructure(
            "no sulfonamide (-SO2NH2) group found; this module only "
            "handles sulfonamides"
        )
    if len(sulfur_atoms) > 1:
        raise UnsupportedStructure(
            "more than one sulfonamide group is out of scope for this "
            "module"
        )
    sulfonamide_atom_idxs = set()
    for s in sulfur_atoms:
        sulfonamide_atom_idxs.add(s.GetIdx())
        sulfonamide_atom_idxs.update(
            n.GetIdx() for n in s.GetNeighbors() if n.GetAtomicNum() in (7, 8)
        )

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
        elif atom.GetIdx() not in sulfonamide_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a sulfonamide group (P-65.3.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular an N-substituted sulfonamide, or a "
                "coexisting carboxylic acid or other characteristic group, "
                "needs handling not yet implemented here"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (sulfur,) = sulfur_atoms
    (carbon,) = (n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 6)
    (nitrogen,) = (n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 7)
    n_alkyl_carbons = tuple(n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != sulfur.GetIdx())
    return sulfur.GetIdx(), carbon.GetIdx(), n_alkyl_carbons


def _reject_enesulfonamide_carbon(graph, so2nh2_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if so2nh2_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a sulfonamide on a carbon that is also part of a C=C/C#C bond "
            "is out of scope for this module"
        )


def _n_alkyl_info(full_carbon_graph, n_alkyl_carbons, non_single_bond_atoms):
    """Mirrors `_amide.py`'s own N-alkyl handling: each N-substituent must be
    a plain, unbranched, unsubstituted, saturated alkyl chain. Returns
    (n_names, n_substituent_atoms) -- the latter must be excluded from the
    carbon graph before picking the principal chain, since these carbons
    hang off the (excluded) sulfonamide nitrogen rather than off the
    sulfonamide carbon itself."""
    n_names = []
    n_substituent_atoms = set()
    for n_alkyl_c in n_alkyl_carbons:
        n_length = linear_branch(full_carbon_graph, n_alkyl_c, None)
        if n_length is None:
            raise UnsupportedStructure("a branched N-substituent is not supported yet")
        n_atoms = set()
        previous, current = None, n_alkyl_c
        while current is not None:
            n_atoms.add(current)
            neighbors = [n for n in full_carbon_graph[current] if n != previous]
            previous, current = current, (neighbors[0] if neighbors else None)
        if n_atoms & non_single_bond_atoms:
            raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
        n_names.append(alkyl_name(n_length))
        n_substituent_atoms |= n_atoms
    return n_names, n_substituent_atoms


def _n_prefix(n_names):
    if len(n_names) == 2 and n_names[0] == n_names[1]:
        return f"N,N-di{n_names[0]}"
    return "-".join(f"N-{name}" for name in sorted(n_names))


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


def _suffix_body(ene_locants, yne_locants, so2nh2_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append(([so2nh2_locant], "sulfonamide"))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, so2nh2_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "sulfonamide"

    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole -SO2NH2) in total omits the locant, e.g.
        # 'ethanesulfonamide'.
        return alkane_name(2) + "sulfonamide"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, so2nh2_locant)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, so2nh2_locant, ene_locants, yne_locants, substituents):
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, so2nh2_locant, ene_locants, yne_locants, grouped)
    return (
        (
            so2nh2_locant,
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


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _substituents_for_ring(graph, ring_order, halogens, excluded):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, so2nh2_locant, grouped):
    stem = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if total_subs == 0:
        # P-14.3.3: the sole substituent on an otherwise unsubstituted ring
        # has no locant to distinguish, e.g. 'cyclohexanesulfonamide'.
        return stem + "sulfonamide"

    prefix = format_substituent_prefixes(grouped)
    return f"{prefix}{stem}-{so2nh2_locant}-sulfonamide"


def _ring_candidate_key(ring_size, so2nh2_locant, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _ring_name_from_substituents(ring_size, so2nh2_locant, grouped)
    return so2nh2_locant, locant_set, citation_locants, name


def _name_cyclic_sulfonamide(mol, sulfur_idx, so2nh2_carbon, n_names, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must lie on
    the ring itself (P-92: a stereocenter on a substituent branch is out
    of scope, mirroring `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid`),
    and the winning ring numbering's own locants for those atoms are used
    to format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
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
            so2nh2_locant = position_of[so2nh2_carbon]
            substituents = _substituents_for_ring(graph, candidate, halogens, excluded)
            key = _ring_candidate_key(ring_size, so2nh2_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if n_names:
        best_name = f"{_n_prefix(n_names)}{best_name}"

    if stereo is not None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def name_sulfonamide(mol) -> str:
    sulfur_idx, so2nh2_carbon, n_alkyl_carbons = _validate_and_collect_sulfonamides(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and sulfur_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 2:
        # The two S=O double bonds are always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enesulfonamide_carbon(graph, so2nh2_carbon, bonds)
    non_single_bond_atoms = {a for a, b, _ in all_non_single} | {b for a, b, _ in all_non_single}
    n_names, n_substituent_atoms = _n_alkyl_info(carbon_adjacency(mol), n_alkyl_carbons, non_single_bond_atoms)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        raise UnsupportedStructure(
            "polycyclic/spiro sulfonamides are not supported yet (this "
            "module only handles acyclic chains and a single saturated "
            "ring)"
        )
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        ring_atoms = set(ring_info.AtomRings()[0])
        if so2nh2_carbon not in ring_atoms:
            raise UnsupportedStructure(
                "a sulfonamide on a substituent branch chain rather than "
                "the ring itself is not supported yet"
            )
        return _name_cyclic_sulfonamide(mol, sulfur_idx, so2nh2_carbon, n_names, stereo)

    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    full_carbon_graph = carbon_adjacency(mol)
    carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in n_substituent_atoms}
    chains = _longest_chains(carbon_graph)
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if so2nh2_carbon not in chain:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        if stereo is not None and any(atom not in chain for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            so2nh2_carbon in c and (not bonds or _bond_locants(c, bonds) is not None) for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the sulfonamide-bearing carbon (and/or a multiple bond) does "
            "not lie on a single longest carbon chain; a shorter principal "
            "chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2nh2_locant = position_of[so2nh2_carbon]
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, so2nh2_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if n_names:
        best_name = f"{_n_prefix(n_names)}{best_name}"

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
