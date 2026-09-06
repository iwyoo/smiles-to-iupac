"""Naming of tellurols (the '-tellurol' suffix, -TeH) on acyclic saturated
or unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-63.1.1, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'tellurol' is the next chalcogen analogue after 'selenol' (`_selenol.py`)
  in the same family as 'thiol' (`_thiol.py`) — 'ol'/'thiol'/'selenol'/
  'tellurol' are all the same substitutive-suffix construction, just with
  O/S/Se/Te respectively. This module mirrors `_selenol.py` exactly,
  swapping the chalcogen atomic number (52 instead of 34) and suffix word
  ('tellurol' instead of 'selenol') and otherwise reusing the same
  `_common.py` chain-search/locant helpers directly. Confirmed via
  PubChem PUG REST: CID 356643 (`C[TeH]`) -> "methanetellurol", CID
  71407255 (`CC[TeH]`) -> "ethanetellurol" (the same P-14.3.4.2(b)
  two-carbon locant omission as 'ethaneselenol'/'ethanethiol'), CID
  71405781 (`CCC[TeH]`) -> "propane-1-tellurol", CID 101718199
  (`CC(C)C[TeH]`) -> "2-methylpropane-1-tellurol" (branched), CID
  101841305 (`C=CC[TeH]`) -> "prop-2-ene-1-tellurol" (unsaturated chain).
- P-35.2.1 (Chapter P-3): a halogen substituent on the carbon chain is
  prefix-only and coexists freely with the -TeH suffix, the same
  mechanism already independently verified in `_thiol.py`/`_selenol.py`/
  several other modules -- no PubChem-listed halogenated tellurol was
  found to independently confirm this specific combination (`ClCC[TeH]`
  came back as CID 0, same as `_selenol.py`'s equivalent case), so this
  is a reviewed (eyeballed), not independently verified, extension along
  an otherwise already-confirmed single axis.
- Two or more -TeH groups (ditellurol, tritellurol, ...) lifts the
  group-count cap entirely, the same way `_selenol.py`'s own triselenol+
  generalization did, and for the identical reason: `_suffix_body`/
  `_name_from_substituents` below are already fully generalized over an
  arbitrary-length `te_locants` list, so lifting the cap adds no new code
  path to verify. No ditellurol compound was found registered in PubChem
  at all (`[TeH]CC[TeH]`/`[TeH]CCC[TeH]` both came back as CID 0, an even
  weaker evidence bar than `_selenol.py`'s triselenol case, which at
  least found a registered but IUPACName-less structure) -- this is a
  reviewed, not independently structure-verified, generalization along an
  already-confirmed axis, matching `_thiol.py`'s and `_selenol.py`'s own
  precedent.
- A single, otherwise-unsubstituted saturated monocyclic ring also works,
  mirroring `_thiol.py`'s/`_selenol.py`'s own monocyclic support:
  confirmed via PubChem structure match, `C1CCCCC1[TeH]` ->
  "cyclohexanetellurol". The ring machinery below is ported directly from
  `_selenol.py`'s own ring code (itself ported from `_thiol.py`), so a
  multi-tellurol ring is supported too on the same trust-the-shared-
  machinery basis as `_selenol.py`'s own multi-group ring support.

- P-91.3/P-92: a
  molecule with one or more *specified* tetrahedral stereocenters -- every
  one on the principal chain/ring itself, no unspecified one alongside
  them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
  `_thiol.py`/`_selenol.py` (chain and ring both). Like -SH/-SeH, a
  tellurol's -TeH tellurium is monovalent (bonded only to its one carbon
  and one H) and confirmed via RDKit's `Chem.FindPotentialStereo` to
  never itself be a potential stereocenter, so this support is
  unconditional.

Scope, deliberately narrow (mirrors `_selenol.py`'s own group-count- and
ring-generalized scope): one or more -TeH groups on an acyclic chain, or
on a single saturated carbon ring, with no other heteroatom (in
particular no -OH, -SH, -SeH, or amine nitrogen) anywhere in the
molecule.

P-31.1.3: a monocyclic ring bearing a tellurol and exactly one C=C ring
double bond -- e.g. 'cyclohex-2-ene-1-tellurol' -- mirrors
`_thiol.py`/`_selenol.py`'s identical extension; the exact ring
structures aren't PubChem-registered (same sparse-coverage reason as the
plain 'cyclohexanetellurol' case above), but the acyclic ene+tellurol
combination is already PubChem-confirmed ('prop-2-ene-1-tellurol', CID
101841305). Deliberately narrow: a ring triple bond, and any other
substituent alongside the ring double bond, are both still explicitly
rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): polycyclic/spiro
rings, unsaturation reaching outside the ring or a ring triple bond, a
-TeH on a substituent branch off an otherwise-unsubstituted ring, a
telluride (-Te- ether-analogue) or any other tellurium-oxidation-state
group, a thiol/selenol or other chalcogen atom, and any oxygen or
nitrogen atom at all. `_name_benzenetellurol` names a single -TeH
directly on a benzene ring carbon (with or without other ring
substituents), e.g. 'benzenetellurol' (PubChem CID 5246059), mirroring
`_thiol.py`'s/`_selenol.py`'s identical construction -- the -TeH's own
locant is never cited, unlike the cycloalkane case; two or more -TeH
groups directly on the ring remain out of scope.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ring_chain_attachment,
    ring_cycle,
    specified_stereocenters,
)
from ._numerals import alkane_name, alkyl_name
from ._substituents import alpha_sort_key, branch_atom_locant, format_substituent_prefixes, name_branch

_YNE_BOND_ORDER = 3.0
_TELLURIUM = 52
_ALLOWED_ATOMIC_NUMS = {6, _TELLURIUM, *HALOGEN_PREFIXES}


def has_tellurol_shape(mol) -> bool:
    return any(atom.GetAtomicNum() == _TELLURIUM for atom in mol.GetAtoms())


def _validate_and_collect_tellurols(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_tellurol`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_tellurol`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors `_thiol.py`'s
    `_validate_and_collect_thiols`."""
    tellurols = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a tellurol tellurium (P-63.1.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting -OH, -SH, -SeH, or amine "
                "nitrogen needs Table 3.3 seniority-coexistence handling "
                "not yet implemented for tellurols"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure("aromatic rings are out of scope for this module")
        elif atomic_num == _TELLURIUM:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a tellurium bonded to more than one heavy atom (e.g. "
                    "a telluride) is out of scope; only an isolated "
                    "tellurol (-TeH) is supported (P-63.1.1)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure("a tellurium double-bonded to carbon is not a tellurol")
            if atom.GetTotalNumHs() != 1:
                raise UnsupportedStructure(
                    "a -Te- atom that isn't a simple tellurol (-TeH) is out "
                    "of scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a tellurol must be attached to a carbon atom")
            tellurols.add(atom.GetIdx())
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent hydride to substitute"
        )
    if not tellurols:
        raise UnsupportedStructure(
            "no tellurol (-TeH) group found; this module only handles tellurols"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return tellurols


def _reject_enetellurol_carbons(graph, tellurols, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for te_idx in tellurols:
        (carbon,) = graph[te_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a tellurol on a carbon that is also part of a C=C/C#C bond "
                "is out of scope for this module"
            )


def _suffix_body(ene_locants, yne_locants, te_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(te_locants), multiplied_word(len(te_locants), "tellurol")))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _name_from_substituents(chain_length, te_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants are always '1' and
        # never cited.
        tellurol_word = multiplied_word(len(te_locants), "tellurol")
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + tellurol_word

    if chain_length == 2 and not has_unsaturation and total_subs == 0 and len(te_locants) == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (here, the sole -TeH) in total omits the locant,
        # e.g. 'ethanetellurol'.
        return alkane_name(2) + "tellurol"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, te_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, te_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    te_locant_set = lowest_locant_set(te_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, te_locants, ene_locants, yne_locants, grouped)
    return (
        (
            te_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _te_locants(position_of, tellurols, graph):
    locants = []
    for te in tellurols:
        (carbon,) = graph[te]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, tellurols):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in tellurols]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _substituents_for_ring(graph, ring_order, halogens, tellurols):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in tellurols]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, te_locants, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    stem = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if not has_unsaturation:
        tellurol_word = multiplied_word(len(te_locants), "tellurol")
        if total_subs == 0 and len(te_locants) == 1:
            # P-14.3.3: the sole substituent on an otherwise unsubstituted
            # ring has no locant to distinguish, e.g. 'cyclohexanetellurol'.
            return stem + tellurol_word
        prefix = format_substituent_prefixes(grouped)
        loc_str = ",".join(str(loc) for loc in sorted(te_locants))
        return f"{prefix}{stem}-{loc_str}-{tellurol_word}"

    # A competing ring double/triple bond (P-31.1.3) means the tellurol's
    # locant is never omittable even when it's the sole substituent --
    # mirrors `_thiol.py`/`_selenol.py`'s identical treatment.
    unsaturated_stem = stem[:-3]
    needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    prefix = format_substituent_prefixes(grouped)
    body = _suffix_body(ene_locants, yne_locants, te_locants)
    return prefix + unsaturated_stem + ("a" if needs_stem_a else "") + "-" + body


def _ring_candidate_key(ring_size, te_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    te_locant_set = lowest_locant_set(te_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, te_locants, ene_locants, yne_locants, grouped)
    return te_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_bond_locant(position_of, bond_atoms, ring_size):
    pa, pb = position_of[bond_atoms[0]], position_of[bond_atoms[1]]
    return ring_size if {pa, pb} == {1, ring_size} else min(pa, pb)


def _ring_bond_locants(position_of, bonds, ring_size):
    """(ene_locants, yne_locants), both sorted, for every ring C=C/C#C bond
    under this ring numbering -- mirrors `_thiol.py`/`_selenol.py`'s
    identical helper."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _ring_bond_locant(position_of, (a, b), ring_size)
        (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return sorted(ene), sorted(yne)


def _ring_branch_stereo_display(graph, ring_order, tellurols, stereo, halogens):
    """Mirrors `_thiol.py`/`_selenol.py`'s identical helper (itself
    mirroring `_aromatic.py`'s `_stereo_display`): if the ring carries
    exactly one specified stereocenter and that stereocenter sits off the
    ring on the ring's own sole substituent branch (P-92), return that
    branch's ring-attachment atom plus its bracketed
    "[(<locant><R/S>)-<name>]" display (P-91.3). Returns None (the caller
    keeps its existing outright rejection) for more than one stereocenter,
    or the ring having more or fewer than one substituent in total."""
    if len(stereo) != 1:
        return None
    stereo_atom, r_or_s = stereo[0]
    ring_set = set(ring_order)
    branch_attachments = [
        (ring_atom, neighbor)
        for ring_atom in ring_order
        for neighbor in graph[ring_atom]
        if neighbor not in ring_set and neighbor not in tellurols
    ]
    if len(branch_attachments) != 1:
        return None
    ring_atom, branch_root = branch_attachments[0]
    branch_name, branch_compound = name_branch(graph, branch_root, ring_atom, halogens)
    site_locant = branch_atom_locant(graph, branch_root, ring_atom, stereo_atom, halogens)
    descriptor = f"({site_locant}{r_or_s})-{branch_name}"
    display = f"[{descriptor}]" if branch_compound else f"({descriptor})"
    return ring_atom, display


def _name_cyclic_tellurol(mol, tellurols, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_thiol.py`'s `_name_cyclic_thiol`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_thiol.py`/`_selenol.py`'s
    own case): exactly one stereocenter on the ring's sole substituent
    branch instead embeds a bracketed descriptor into that substituent's
    own name, in place of the usual ring-locant prefix."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, tellurols, stereo, halogens)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(_substituents_for_ring(graph, ring_order, halogens, tellurols).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "tellurol is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            te_locants = _te_locants(position_of, tellurols, graph)
            substituents = _substituents_for_ring(graph, candidate, halogens, tellurols)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = _ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, te_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _benzenetellurol_name_from_substituents(te_locants, grouped):
    # Unlike the cycloalkane case, the mancude ring's own numbering is
    # always free to start at the -TeH carbon (P-14.3.3-style), so its
    # locant is never cited even when other substituents need theirs,
    # mirroring `_thiol.py`'s/`_selenol.py`'s identical 'benzenethiol'/
    # 'benzeneselenol' treatment.
    tellurol_word = multiplied_word(len(te_locants), "tellurol")
    if not grouped:
        return "benzene" + tellurol_word
    return f"{format_substituent_prefixes(grouped)}benzene{tellurol_word}"


def _benzenetellurol_candidate_key(te_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    te_locant_set = lowest_locant_set(te_locants)
    name = _benzenetellurol_name_from_substituents(te_locants, grouped)
    return te_locant_set, locant_set, citation_locants, name


def _name_benzenetellurol(mol, ring_atoms):
    """P-63.1.1: -TeH attached directly to a benzene ring carbon -- e.g.
    'benzenetellurol' (PubChem CID 5246059). Mirrors `_thiol.py`'s/
    `_selenol.py`'s `_name_benzenethiol`/`_name_benzeneselenol` exactly
    (tellurium in place of sulfur/selenium), with the retained name
    'benzene' as stem in place of 'cyclo' + alkane_name. Only a single
    -TeH directly on the ring is verified here (two or more direct ring
    tellurols remain out of scope)."""
    tellurols = _validate_and_collect_tellurols(mol, aromatic_ring_atoms=ring_atoms)
    if len(tellurols) != 1:
        raise UnsupportedStructure(
            "more than one tellurol directly on the benzene ring is not "
            "supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside benzenetellurol is not "
            "supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            te_locants = _te_locants(position_of, tellurols, graph)
            substituents = _substituents_for_ring(graph, candidate, halogens, tellurols)
            key = _benzenetellurol_candidate_key(te_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_tellurol(mol, ring_atoms):
    """Name a tellurol whose -TeH lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 3-phenylpropane-1-tellurol. The ring is cited as
    a 'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_thiol.py`'s `_name_phenyl_chain_thiol`. Narrower than the acyclic
    path above: exactly one -TeH, no chain unsaturation, and no specified
    stereocenter."""
    tellurols = _validate_and_collect_tellurols(mol, aromatic_ring_atoms=ring_atoms)
    if len(tellurols) != 1:
        raise UnsupportedStructure(
            "more than one tellurol alongside a benzene-ring substituent "
            "is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "tellurol chain is not supported yet"
        )
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "tellurol chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain tellurol is not supported yet"
        )
    ring_atom, chain_root = attachment
    if chain_root in tellurols:
        raise UnsupportedStructure(
            "a tellurol directly on the benzene ring (tellurophenol-type) "
            "uses a separate construction, out of scope for this "
            "chain-parent module"
        )
    anchor_tellurium = next(iter(tellurols))
    (anchor_carbon,) = graph[anchor_tellurium]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, tellurols)
    chain_set = set(chain)
    for t in tellurols:
        (carbon,) = graph[t]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a tellurol outside the single unbranched chain hanging off "
                "the benzene ring is not supported yet"
            )
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    halogens = halogen_substituents(mol)
    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        te_locants = _te_locants(position_of, tellurols, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, te_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_tellurol(mol, tellurols):
    """Name one or more -TeH groups lying entirely on a single branched
    chain hanging off one atom of an otherwise-plain saturated monocyclic
    ring (the ring itself bears no tellurol) -- e.g.
    cyclohexylmethanetellurol. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_tellurol` above and
    `_thiol.py`'s `_name_ring_substituent_chain_thiol`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, tellurols)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    anchor_tellurium = next(iter(tellurols))
    (anchor_carbon,) = graph[anchor_tellurium]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, tellurols)
    chain_set = set(chain)
    for t in tellurols:
        (carbon,) = graph[t]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a tellurol outside the single branched chain hanging off "
                "the ring is not supported yet"
            )
    branches_by_atom = {
        chain[position - 1]: [r for r in roots if r != ring_atom]
        for position, roots in branches.items()
    }
    branches_by_atom = {atom: roots for atom, roots in branches_by_atom.items() if roots}

    ring_name = "cyclo" + alkyl_name(len(ring_atoms))
    chain_length = len(chain)

    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        te_locants = _te_locants(position_of, tellurols, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, te_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_tellurol(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            ring_tellurols = _validate_and_collect_tellurols(mol, aromatic_ring_atoms=ring_atoms)
            if len(ring_tellurols) == 1:
                (only_te,) = ring_tellurols
                (only_te_carbon,) = adjacency(mol)[only_te]
                if only_te_carbon in ring_atoms:
                    return _name_benzenetellurol(mol, ring_atoms)
            return _name_phenyl_chain_tellurol(mol, ring_atoms)
    tellurols = _validate_and_collect_tellurols(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, _YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enetellurol_carbons(graph, tellurols, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        raise UnsupportedStructure(
            "polycyclic/spiro tellurols are not supported yet (this module "
            "only handles acyclic chains and a single saturated ring)"
        )
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic tellurol "
                "is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == _YNE_BOND_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a tellurol is "
                "not supported yet -- only a ring double bond is in scope "
                "for this first pass (see P-31.1.3)"
            )
        ring_tellurols = {t for t in tellurols if next(iter(graph[t])) in ring_atoms}
        if not bonds and not ring_tellurols:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than "
                    "the ring itself is not supported yet (see P-92)"
                )
            return _name_ring_substituent_chain_tellurol(mol, tellurols)
        if ring_tellurols != tellurols:
            raise UnsupportedStructure(
                "a tellurol on a substituent branch chain rather than the "
                "ring itself is not supported yet"
            )
        return _name_cyclic_tellurol(mol, tellurols, stereo, bonds)

    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _te_locants(position_of, tellurols, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _te_locants({a: i + 1 for i, a in enumerate(c)}, tellurols, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every tellurol-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain, "
            "or a -TeH expressed as a 'tellanyl' substituent prefix, is "
            "not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            te_locants = _te_locants(position_of, tellurols, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, tellurols)
            key, name = _candidate_key(chain_length, te_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
