"""Naming of sulfinic acids (the '-sulfinic acid' suffix, -SO2H) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.3.1, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  'sulfinic acid' is the preselected suffix for -S(=O)-OH, one oxidation
  state below sulfonic acid (`_sulfonic_acid.py`) -- the sulfur carries one
  carbon, one double-bonded oxygen, and one hydroxyl oxygen (degree 3, not
  4). Table 3.3's acid cluster ranks sulfinic acid junior to sulfonic acid;
  this module doesn't implement that acid-vs-acid seniority (out of scope,
  see below -- this module only ever handles a single sulfinic acid group
  in isolation).
- Like 'sulfonic acid' (and unlike 'ol'/'thiol'/'amine'), 'sulfinic acid'
  is cited as the parent hydride name followed directly by the two-word
  suffix with no elision -- 'ethane' + 'sulfinic acid' -> 'ethanesulfinic
  acid'.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_sulfonic_acid.py` apply, e.g. 'ethanesulfinic acid'.
- P-44.4.1.8 / P-45.2: the -SO2H locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as every other suffix module here.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -SO2H suffix.
- P-92 stereocenters: unlike
  `_sulfonic_acid.py`'s sulfur (two identical =O, never stereogenic), this
  module's sulfinic sulfur (-R, =O, -OH, a lone pair) is *itself* a
  potential stereocenter in essentially every real -SO2H molecule --
  confirmed via RDKit's `Chem.FindPotentialStereo` on e.g. plain
  `CC(C)S(=O)O` (no chain stereocenter at all), which still flags the
  sulfur atom. That means any input with a specified *chain* stereocenter
  almost always has this second, unspecified sulfur stereocenter riding
  along, and `_common.specified_stereocenters` correctly rejects that
  combination as partially specified (P-92) rather than silently ignoring
  either one. A specified sulfur configuration is out of scope too, since
  this project has no established locant/prefix convention for a
  heteroatom-centered (rather than carbon-centered) stereodescriptor, and
  PubChem itself doesn't distinguish the two sulfur configurations of a
  test case (`CCC[S@](=O)O`/`CCC[S@@](=O)O`, both CID 643586, same
  unstereo name) -- so this module only ever explicitly rejects a
  specified stereocenter (chain carbon or sulfur alike) rather than
  attempting to cite one; see the module below for the R/S support this
  project *does* provide (`_carboxylic_acid.py`/`_aldehyde.py`/
  `_ketone.py`/`_sulfonic_acid.py`, none of which have this complication).

Scope, deliberately narrow (mirrors `_sulfonic_acid.py`'s own first pass):
a single -SO2H on an acyclic chain or on a single saturated carbon ring
(monocyclic, mirroring `_sulfonic_acid.py`'s own monocyclic support),
with no other heteroatom anywhere in the molecule except the sulfinic
acid group's own two oxygens -- acid-vs-acid seniority and any other
Table 3.3 seniority coexistence is future work, tracked under
`multi-carbonyl-seniority.md`.

P-31.1.3: a monocyclic ring bearing a sulfinic acid and exactly one C=C
ring double bond -- e.g. 'cyclohex-2-ene-1-sulfinic acid',
'cyclohex-3-ene-1-sulfinic acid', both confirmed via PubChem PUG REST.
Mirrors `_sulfonic_acid.py`'s identical extension; deliberately narrow: a
ring triple bond, and any other substituent alongside the ring double
bond, are both still explicitly rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): polycyclic/spiro
rings, unsaturation reaching outside the ring or a ring triple bond, an
-SO2H on a substituent branch off an otherwise-unsubstituted *saturated*
ring, two or more -SO2H groups, and a sulfinic acid on a carbon that is
also part of a C=C/C#C bond. Two *aromatic*-ring cases:
`_name_phenyl_chain_sulfinic_acid` names a -SO2H chain hanging off a
single plain, unsubstituted benzene ring (e.g.
'3-phenylpropane-1-sulfinic acid'), mirroring `_sulfonic_acid.py`'s
identical benzene-ring-substituent path -- narrower than the acyclic
path: no chain unsaturation. `_name_benzenesulfinic_acid` names -SO2H
directly on a benzene ring carbon (with or without other ring
substituents), e.g. 'benzenesulfinic acid' (PubChem CID 12057),
'2-methylbenzenesulfinic acid' (CID 12661295), mirroring
`_sulfonic_acid.py`'s `_name_benzenesulfonic_acid` with the retained
name 'benzene' as stem -- the -SO2H's own locant is never cited here,
unlike the cycloalkane case.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    halogen_substituents,
    is_plain_benzene_ring,
    lowest_locant_set,
    non_single_bonds,
    ordered_chain,
    path_between,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    ring_cycle,
    specified_stereocenters,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import (
    alpha_sort_key,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
)

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _sulfinic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfinic acid group: bonded to exactly one
    carbon, one double-bonded (terminal) oxygen, and one single-bonded
    hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 2:
            continue
        double_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0]
        hydroxyl_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(double_os) != 1 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_sulfinic_acid_shape(mol) -> bool:
    return bool(_sulfinic_sulfur_atoms(mol))


def _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_sulfinic_acid`'s benzene-
    ring-substituent path (see `_name_phenyl_chain_sulfinic_acid`) can
    reuse this same validation for the rest of the molecule. Empty by
    default, so every other caller's behavior is unchanged."""
    sulfur_atoms = _sulfinic_sulfur_atoms(mol)
    if not sulfur_atoms:
        raise UnsupportedStructure(
            "no sulfinic acid (-SO2H) group found; this module only "
            "handles sulfinic acids"
        )
    if len(sulfur_atoms) > 1:
        raise UnsupportedStructure(
            "more than one sulfinic acid group is out of scope for this "
            "module"
        )
    sulfinic_atom_idxs = set()
    for s in sulfur_atoms:
        sulfinic_atom_idxs.add(s.GetIdx())
        sulfinic_atom_idxs.update(n.GetIdx() for n in s.GetNeighbors() if n.GetAtomicNum() == 8)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
        elif atom.GetIdx() not in sulfinic_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a sulfinic acid group (P-65.3.1) "
                "and halogen substituents (P-35.2.1) are not supported "
                "yet -- in particular a coexisting carboxylic/sulfonic "
                "acid or other characteristic group needs acid-vs-acid "
                "Table 3.3 seniority handling not yet implemented here"
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
    return sulfur.GetIdx(), carbon.GetIdx()


def _reject_enesulfinic_carbon(graph, so2h_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if so2h_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a sulfinic acid on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
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


def _suffix_body(ene_locants, yne_locants, so2h_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append(([so2h_locant], "sulfinic acid"))

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


def _name_from_substituents(chain_length, so2h_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "sulfinic acid"

    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole -SO2H) in total omits the locant, e.g.
        # 'ethanesulfinic acid'.
        return alkane_name(2) + "sulfinic acid"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, so2h_locant)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, so2h_locant, ene_locants, yne_locants, substituents):
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
    name = _name_from_substituents(chain_length, so2h_locant, ene_locants, yne_locants, grouped)
    return (
        (
            so2h_locant,
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


def _ring_name_from_substituents(ring_size, so2h_locant, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    stem = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if not has_unsaturation:
        if total_subs == 0:
            # P-14.3.3: the sole substituent on an otherwise unsubstituted
            # ring has no locant to distinguish, e.g.
            # 'cyclohexanesulfinic acid'.
            return stem + "sulfinic acid"
        prefix = format_substituent_prefixes(grouped)
        return f"{prefix}{stem}-{so2h_locant}-sulfinic acid"

    # A competing ring double/triple bond (P-31.1.3) means the sulfinic
    # acid's locant is never omittable even when it's the sole
    # substituent -- mirrors `_sulfonic_acid.py`'s identical treatment.
    unsaturated_stem = stem[:-3]
    needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    prefix = format_substituent_prefixes(grouped)
    body = _suffix_body(ene_locants, yne_locants, so2h_locant)
    return prefix + unsaturated_stem + ("a" if needs_stem_a else "") + "-" + body


def _ring_candidate_key(ring_size, so2h_locant, ene_locants, yne_locants, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, so2h_locant, ene_locants, yne_locants, grouped)
    return so2h_locant, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_bond_locant(position_of, bond_atoms, ring_size):
    pa, pb = position_of[bond_atoms[0]], position_of[bond_atoms[1]]
    return ring_size if {pa, pb} == {1, ring_size} else min(pa, pb)


def _ring_bond_locants(position_of, bonds, ring_size):
    """(ene_locants, yne_locants), both sorted, for every ring C=C/C#C bond
    under this ring numbering -- mirrors `_sulfonic_acid.py`'s identical
    helper."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _ring_bond_locant(position_of, (a, b), ring_size)
        (ene if order == _ENE_ORDER else yne).append(locant)
    return sorted(ene), sorted(yne)


def _name_cyclic_sulfinic_acid(mol, sulfur_idx, so2h_carbon, bonds=()):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    if bonds and any(_substituents_for_ring(graph, ring_order, halogens, excluded).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "sulfinic acid is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2h_locant = position_of[so2h_carbon]
            substituents = _substituents_for_ring(graph, candidate, halogens, excluded)
            ene_locants, yne_locants = _ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, so2h_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _benzenesulfinic_acid_name_from_substituents(grouped):
    # Unlike the cycloalkane case, the mancude ring's own numbering is
    # always free to start at the -SO2H carbon (P-14.3.3-style), so its
    # locant is never cited even when other substituents need theirs,
    # e.g. '2-methylbenzenesulfinic acid' (PubChem CID 12661295), not
    # '2-methylbenzene-1-sulfinic acid' -- mirrors `_sulfonic_acid.py`'s
    # identical 'benzenesulfonic acid' treatment.
    if not grouped:
        return "benzenesulfinic acid"
    return f"{format_substituent_prefixes(grouped)}benzenesulfinic acid"


def _benzenesulfinic_acid_candidate_key(so2h_locant, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _benzenesulfinic_acid_name_from_substituents(grouped)
    return so2h_locant, locant_set, citation_locants, name


def _name_benzenesulfinic_acid(mol, ring_atoms):
    """P-65.3.1: -SO2H attached directly to a benzene ring carbon -- e.g.
    'benzenesulfinic acid' (PubChem CID 12057), '2-methylbenzenesulfinic
    acid' (CID 12661295). Mirrors `_sulfonic_acid.py`'s
    `_name_benzenesulfonic_acid` exactly, with the retained name
    'benzene' as stem in place of 'cyclo' + alkane_name; an aromatic ring
    has no ene/yne locants of its own."""
    sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol) is not None:
        # See module docstring: the sulfinic sulfur is itself a potential
        # stereocenter in virtually every real -SO2H molecule, and this
        # project has no established way to cite one -- reject
        # unconditionally, same as the acyclic/chain paths.
        raise UnsupportedStructure(
            "a specified stereocenter (ring carbon or the sulfinic "
            "sulfur itself) is not supported yet for sulfinic acids (see "
            "P-92, module docstring)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2h_locant = position_of[so2h_carbon]
            substituents = _substituents_for_ring(graph, candidate, halogens, excluded)
            key = _benzenesulfinic_acid_candidate_key(so2h_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_sulfinic_acid(mol, ring_atoms):
    """Name a sulfinic acid whose -SO2H lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropane-1-sulfinic acid.
    The ring is cited as a 'phenyl' substituent prefix (via
    `name_branch`'s aromatic-ring recognition) on the chain, which is the
    parent hydride, mirroring `_sulfonic_acid.py`'s
    `_name_phenyl_chain_sulfonic_acid`. Narrower than the acyclic path
    above: no chain unsaturation -- a separate follow-up (see
    tasks/phenyl-substituent-on-sulfinic-acid-chain.md's scope note)."""
    sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol) is not None:
        # See module docstring: the sulfinic sulfur is itself a potential
        # stereocenter in virtually every real -SO2H molecule, and this
        # project has no established way to cite one -- reject
        # unconditionally, same as the acyclic path below.
        raise UnsupportedStructure(
            "a specified stereocenter (chain carbon or the sulfinic "
            "sulfur itself) is not supported yet for sulfinic acids (see "
            "P-92, module docstring)"
        )
    excluded = {sulfur_idx}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if sulfur_idx not in (b[0], b[1]) and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "sulfinic acid chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain sulfinic acid is "
            "not supported yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, excluded)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside a "
            "sulfinic acid is not supported yet"
        )
    if so2h_carbon not in chain:
        raise UnsupportedStructure(
            "the sulfinic acid carbon must lie on the chain hanging off "
            "the benzene ring for this benzene-substituent path"
        )

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        so2h_locant = position_of[so2h_carbon]
        substituents = {
            position_of[chain_root]: [name_branch(graph, ring_atom, chain_root, halogens, ring_atoms)]
        }
        key, name = _candidate_key(chain_length, so2h_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_sulfinic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=ring_atoms)
            if so2h_carbon in ring_atoms:
                return _name_benzenesulfinic_acid(mol, ring_atoms)
            return _name_phenyl_chain_sulfinic_acid(mol, ring_atoms)
    sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol)
    if specified_stereocenters(mol) is not None:
        # Unlike `_sulfonic_acid.py`'s sulfur, this module's sulfinic
        # sulfur is itself a potential stereocenter in virtually every real
        # -SO2H molecule (module docstring), and this project has no
        # established way to cite a heteroatom-centered stereodescriptor --
        # explicitly reject rather than silently drop the marker (P-92).
        raise UnsupportedStructure(
            "a specified stereocenter (chain carbon or the sulfinic sulfur "
            "itself) is not supported yet for sulfinic acids (see P-92, "
            "module docstring)"
        )
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and sulfur_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 1:
        # The one S=O double bond is always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enesulfinic_carbon(graph, so2h_carbon, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        raise UnsupportedStructure(
            "polycyclic/spiro sulfinic acids are not supported yet (this "
            "module only handles acyclic chains and a single saturated "
            "ring)"
        )
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic sulfinic "
                "acid is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == _YNE_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a sulfinic "
                "acid is not supported yet -- only a ring double bond is "
                "in scope for this first pass (see P-31.1.3)"
            )
        if so2h_carbon not in ring_atoms:
            raise UnsupportedStructure(
                "a sulfinic acid on a substituent branch chain rather "
                "than the ring itself is not supported yet"
            )
        return _name_cyclic_sulfinic_acid(mol, sulfur_idx, so2h_carbon, bonds)

    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if so2h_carbon not in chain:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the sulfinic-acid-bearing carbon (and/or a multiple bond) "
            "does not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2h_locant = position_of[so2h_carbon]
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, so2h_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
