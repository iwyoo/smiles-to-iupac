"""Naming of selenols (the '-selenol' suffix, -SeH) on acyclic saturated or
unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-63.1.1, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'selenol' is the next chalcogen analogue in the same family as 'thiol'
  (`_thiol.py`) — 'ol'/'thiol'/'selenol'/'tellurol' are all the same
  substitutive-suffix construction, just with O/S/Se/Te respectively.
  This module mirrors `_thiol.py`'s own first-pass scope (single -SH,
  before that module's later dithiol/monocyclic extensions) exactly,
  swapping the chalcogen atomic number (34 instead of 16) and suffix word
  ('selenol' instead of 'thiol') and otherwise reusing the same
  chain-search/locant machinery from `_common.py` directly (this module is
  new enough to use the shared `longest_chains`/`bond_locant`/
  `bond_locants`/`group_substituents`/`multiplied_word` helpers rather
  than adding yet another private copy of them -- `_thiol.py` itself still
  has its own pre-`_common.py` copies, a known, separate cleanup debt not
  addressed here). Confirmed via PubChem PUG REST: CID 440764 (`C[SeH]`)
  -> "methaneselenol", CID 5252527 (`CC[SeH]`) -> "ethaneselenol" (the
  same P-14.3.4.2(b) two-carbon locant omission as 'ethanethiol'), CID
  71373846 (`CCC[SeH]`) -> "propane-1-selenol", CID 157368643
  (`CC(C)C[SeH]`) -> "2-methylpropane-1-selenol" (branched), CID 15821407
  (`C=CC[SeH]`) -> "prop-2-ene-1-selenol" (unsaturated chain).
- P-35.2.1 (Chapter P-3): a halogen substituent on the carbon chain is
  prefix-only and coexists freely with the -SeH suffix, the same
  mechanism already independently verified in `_thiol.py`/`_alcohol.py`/
  several other modules -- no PubChem-listed halogenated selenol was
  found to independently confirm this specific combination (every
  candidate SMILES tried came back as CID 0), so this is a reviewed
  (eyeballed), not independently verified, extension along an otherwise
  already-confirmed single axis.
- Two -SeH groups (a diselenol) mirrors `_thiol.py`'s dithiol extension --
  the locant/suffix machinery already generalizes over a list of
  selenol locants, so only the single-group guard needed lifting.
  Confirmed via PubChem PUG REST: `[SeH]CC[SeH]` ->
  "ethane-1,2-diselenol", `[SeH]CCC[SeH]` -> "propane-1,3-diselenol".
- Three or more -SeH groups (triselenol, tetraselenol, ...) lifts the
  group-count cap entirely, the same way `_thiol.py` itself generalizes
  to "dithiol/trithiol/..." without citing a specific 3+-group PubChem
  worked example of its own -- `_suffix_body`/`_name_from_substituents`
  below are already fully generalized over an arbitrary-length
  `se_locants` list (the same shared mechanism `_alcohol.py`'s own
  confirmed polyol support, e.g. glycerol, uses), so lifting the cap adds
  no new code path to verify. No triselenol compound (three -SeH groups
  on a single chain) was found registered in PubChem to independently
  confirm this specific chalcogen (unlike `_thiol.py`'s own trithiol
  case, which is in the same boat) -- this is a reviewed, not
  independently structure-verified, generalization along an
  already-confirmed axis, matching `_thiol.py`'s own precedent exactly.
- A single, otherwise-unsubstituted saturated monocyclic ring also works,
  mirroring `_thiol.py`'s own monocyclic support: confirmed via PubChem
  structure match, `C1CCCCC1[SeH]` -> "cyclohexaneselenol",
  `C1CCCC1[SeH]` -> "cyclopentaneselenol". The ring machinery below
  (`_ring_name_from_substituents`/`_ring_candidate_key`) is ported
  directly from `_thiol.py`'s own already-generalized-over-multiple-
  groups ring code, so a multi-selenol ring is supported too, the same
  way `_thiol.py` itself doesn't cite a specific multi-group ring PubChem
  example either -- trusting the shared, already-exercised locant
  machinery rather than re-deriving it.

- P-91.3/P-92: a
  molecule with one or more *specified* tetrahedral stereocenters -- every
  one on the principal chain/ring itself, no unspecified one alongside
  them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
  `_thiol.py` (chain and ring both). Like a thiol's -SH sulfur, a
  selenol's -SeH selenium is monovalent (bonded only to its one carbon
  and one H) and confirmed via RDKit's `Chem.FindPotentialStereo` to
  never itself be a potential stereocenter, so this support is
  unconditional.

Scope, deliberately narrow (mirrors `_thiol.py`'s own group-count- and
ring-generalized scope): one or more -SeH groups on an acyclic chain, or
on a single saturated carbon ring, with no other heteroatom (in
particular no -OH, -SH, or amine nitrogen) anywhere in the molecule.

P-31.1.3: a monocyclic ring bearing a selenol and exactly one C=C ring
double bond -- e.g. 'cyclohex-2-ene-1-selenol' -- mirrors `_thiol.py`'s
identical extension (PR #350); the exact ring structures aren't
PubChem-registered (selenol ring coverage there is sparse, same reason
the plain 'cyclohexaneselenol' case above is confirmed only by structure
match rather than name), but the acyclic ene+selenol combination and its
elision rule are already PubChem-confirmed ('prop-2-ene-1-selenol', CID
15821407). Deliberately narrow: a ring triple bond, and any other
substituent alongside the ring double bond, are both still explicitly
rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): polycyclic/spiro
rings, unsaturation reaching outside the ring or a ring triple bond, an
-SeH on a substituent branch off an otherwise-unsubstituted ring, a
selenide (-Se- ether-analogue) or any other selenium-oxidation-state
group, a tellurol or other chalcogen atom, and any oxygen or nitrogen
atom at all. `_name_benzeneselenol` names a single -SeH directly on a
benzene ring carbon (with or without other ring substituents), e.g.
'benzeneselenol' (PubChem CID 69530), mirroring `_thiol.py`'s
`_name_benzenethiol` exactly -- the -SeH's own locant is never cited,
unlike the cycloalkane case; two or more -SeH groups directly on the
ring remain out of scope.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locants,
    carbon_adjacency,
    elides_before,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ordered_chain,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_cycle,
    specified_stereocenters,
)
from ._numerals import alkane_name, alkyl_name
from ._substituents import alpha_sort_key, branch_atom_locant, format_substituent_prefixes, name_branch

_YNE_BOND_ORDER = 3.0
_SELENIUM = 34
_ALLOWED_ATOMIC_NUMS = {6, _SELENIUM, *HALOGEN_PREFIXES}


def has_selenol_shape(mol) -> bool:
    return any(atom.GetAtomicNum() == _SELENIUM for atom in mol.GetAtoms())


def _validate_and_collect_selenols(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_selenol`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_selenol`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors `_thiol.py`'s
    `_validate_and_collect_thiols`."""
    selenols = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a selenol selenium (P-63.1.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting -OH, -SH, or amine nitrogen "
                "needs Table 3.3 seniority-coexistence handling not yet "
                "implemented for selenols"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure("aromatic rings are out of scope for this module")
        elif atomic_num == _SELENIUM:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a selenium bonded to more than one heavy atom (e.g. a "
                    "selenide) is out of scope; only an isolated selenol "
                    "(-SeH) is supported (P-63.1.1)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure("a selenium double-bonded to carbon is not a selenol")
            if atom.GetTotalNumHs() != 1:
                raise UnsupportedStructure(
                    "a -Se- atom that isn't a simple selenol (-SeH) is out "
                    "of scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a selenol must be attached to a carbon atom")
            selenols.add(atom.GetIdx())
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent hydride to substitute"
        )
    if not selenols:
        raise UnsupportedStructure(
            "no selenol (-SeH) group found; this module only handles selenols"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return selenols


def _reject_eneselenol_carbons(graph, selenols, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for se_idx in selenols:
        (carbon,) = graph[se_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a selenol on a carbon that is also part of a C=C/C#C bond "
                "is out of scope for this module"
            )


def _suffix_body(ene_locants, yne_locants, se_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(se_locants), multiplied_word(len(se_locants), "selenol")))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and elides_before(words[i + 1]):
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _name_from_substituents(chain_length, se_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants are always '1' and
        # never cited.
        selenol_word = multiplied_word(len(se_locants), "selenol")
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + selenol_word

    if chain_length == 2 and not has_unsaturation and total_subs == 0 and len(se_locants) == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (here, the sole -SeH) in total omits the locant, e.g.
        # 'ethaneselenol'.
        return alkane_name(2) + "selenol"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, se_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, se_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    se_locant_set = lowest_locant_set(se_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, se_locants, ene_locants, yne_locants, grouped)
    return (
        (
            se_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _se_locants(position_of, selenols, graph):
    locants = []
    for se in selenols:
        (carbon,) = graph[se]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, selenols, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in selenols]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _substituents_for_ring(graph, ring_order, halogens, selenols, mol=None):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in selenols]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, se_locants, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    stem = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if not has_unsaturation:
        selenol_word = multiplied_word(len(se_locants), "selenol")
        if total_subs == 0 and len(se_locants) == 1:
            # P-14.3.3: the sole substituent on an otherwise unsubstituted
            # ring has no locant to distinguish, e.g. 'cyclohexaneselenol'.
            return stem + selenol_word
        prefix = format_substituent_prefixes(grouped)
        loc_str = ",".join(str(loc) for loc in sorted(se_locants))
        return f"{prefix}{stem}-{loc_str}-{selenol_word}"

    # A competing ring double/triple bond (P-31.1.3) means the selenol's
    # locant is never omittable even when it's the sole substituent --
    # mirrors `_thiol.py`'s identical treatment.
    unsaturated_stem = stem[:-3]
    needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    prefix = format_substituent_prefixes(grouped)
    body = _suffix_body(ene_locants, yne_locants, se_locants)
    return prefix + unsaturated_stem + ("a" if needs_stem_a else "") + "-" + body


def _ring_candidate_key(ring_size, se_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    se_locant_set = lowest_locant_set(se_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, se_locants, ene_locants, yne_locants, grouped)
    return se_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, selenols, stereo, halogens, mol=None):
    """Mirrors `_thiol.py`'s identical helper (itself mirroring
    `_aromatic.py`'s `_stereo_display`): if the ring carries exactly one
    specified stereocenter and that stereocenter sits off the ring on the
    ring's own sole substituent branch (P-92), return that branch's
    ring-attachment atom plus its bracketed "[(<locant><R/S>)-<name>]"
    display (P-91.3). Returns None (the caller keeps its existing outright
    rejection) for more than one stereocenter, or the ring having more or
    fewer than one substituent in total."""
    if len(stereo) != 1:
        return None
    stereo_atom, r_or_s = stereo[0]
    ring_set = set(ring_order)
    branch_attachments = [
        (ring_atom, neighbor)
        for ring_atom in ring_order
        for neighbor in graph[ring_atom]
        if neighbor not in ring_set and neighbor not in selenols
    ]
    if len(branch_attachments) != 1:
        return None
    ring_atom, branch_root = branch_attachments[0]
    branch_name, branch_compound = name_branch(graph, branch_root, ring_atom, halogens, mol=mol)
    site_locant = branch_atom_locant(graph, branch_root, ring_atom, stereo_atom, halogens, mol=mol)
    descriptor = f"({site_locant}{r_or_s})-{branch_name}"
    display = f"[{descriptor}]" if branch_compound else f"({descriptor})"
    return ring_atom, display


def _name_cyclic_selenol(mol, selenols, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_thiol.py`'s `_name_cyclic_thiol`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_thiol.py`'s own case):
    exactly one stereocenter on the ring's sole substituent branch instead
    embeds a bracketed descriptor into that substituent's own name, in
    place of the usual ring-locant prefix."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, selenols, stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(_substituents_for_ring(graph, ring_order, halogens, selenols, mol=mol).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "selenol is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            se_locants = _se_locants(position_of, selenols, graph)
            substituents = _substituents_for_ring(graph, candidate, halogens, selenols, mol=mol)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, se_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _benzeneselenol_name_from_substituents(se_locants, grouped):
    # Unlike the cycloalkane case, the mancude ring's own numbering is
    # always free to start at the -SeH carbon (P-14.3.3-style), so its
    # locant is never cited even when other substituents need theirs,
    # mirroring `_thiol.py`'s identical 'benzenethiol' treatment.
    selenol_word = multiplied_word(len(se_locants), "selenol")
    if not grouped:
        return "benzene" + selenol_word
    return f"{format_substituent_prefixes(grouped)}benzene{selenol_word}"


def _benzeneselenol_candidate_key(se_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    se_locant_set = lowest_locant_set(se_locants)
    name = _benzeneselenol_name_from_substituents(se_locants, grouped)
    return se_locant_set, locant_set, citation_locants, name


def _name_benzeneselenol(mol, ring_atoms):
    """P-63.1.1: -SeH attached directly to a benzene ring carbon -- e.g.
    'benzeneselenol' (PubChem CID 69530). Mirrors `_thiol.py`'s
    `_name_benzenethiol` exactly (selenium in place of sulfur), with the
    retained name 'benzene' as stem in place of 'cyclo' + alkane_name.
    Only a single -SeH directly on the ring is verified here (two or more
    direct ring selenols remain out of scope)."""
    selenols = _validate_and_collect_selenols(mol, aromatic_ring_atoms=ring_atoms)
    if len(selenols) != 1:
        raise UnsupportedStructure(
            "more than one selenol directly on the benzene ring is not "
            "supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside benzeneselenol is not "
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
            se_locants = _se_locants(position_of, selenols, graph)
            substituents = _substituents_for_ring(graph, candidate, halogens, selenols, mol=mol)
            key = _benzeneselenol_candidate_key(se_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_selenol(mol, ring_atoms):
    """Name a selenol whose -SeH lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 3-phenylpropane-1-selenol. The ring is cited as
    a 'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_thiol.py`'s `_name_phenyl_chain_thiol`. Narrower than the acyclic
    path above: exactly one -SeH, no chain unsaturation, and no specified
    stereocenter."""
    selenols = _validate_and_collect_selenols(mol, aromatic_ring_atoms=ring_atoms)
    if len(selenols) != 1:
        raise UnsupportedStructure(
            "more than one selenol alongside a benzene-ring substituent is "
            "not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "selenol chain is not supported yet"
        )
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "selenol chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain selenol is not supported yet"
        )
    ring_atom, chain_root = attachment
    if chain_root in selenols:
        raise UnsupportedStructure(
            "a selenol directly on the benzene ring (selenophenol-type) "
            "uses a separate construction, out of scope for this "
            "chain-parent module"
        )
    anchor_selenium = next(iter(selenols))
    (anchor_carbon,) = graph[anchor_selenium]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, selenols, halogens=halogen_substituents(mol))
    chain_set = set(chain)
    for s in selenols:
        (carbon,) = graph[s]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a selenol outside the single unbranched chain hanging off "
                "the benzene ring is not supported yet"
            )
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    halogens = halogen_substituents(mol)
    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        se_locants = _se_locants(position_of, selenols, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, se_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_selenol(mol, selenols):
    """Name one or more -SeH groups lying entirely on a single branched
    chain hanging off one atom of an otherwise-plain saturated monocyclic
    ring (the ring itself bears no selenol) -- e.g.
    cyclohexylmethaneselenol. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_selenol` above and
    `_thiol.py`'s `_name_ring_substituent_chain_thiol`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, selenols)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    anchor_selenium = next(iter(selenols))
    (anchor_carbon,) = graph[anchor_selenium]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, selenols, halogens=halogen_substituents(mol))
    chain_set = set(chain)
    for s in selenols:
        (carbon,) = graph[s]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a selenol outside the single branched chain hanging off "
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
        se_locants = _se_locants(position_of, selenols, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, se_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_with_selenol_chain_selenol(mol, selenols):
    """Name a selenol compound where the ring itself bears at least as
    many -SeH's as a single unbranched chain hanging off exactly one ring
    atom does (P-44.1.1/P-44.1.2.2). The ring is the parent; the chain is
    cited as a '(selanyl...alkyl)' substituent prefix, mirroring
    `_thiol.py`'s `_name_ring_with_thiol_chain_thiol`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, selenols)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, selenols)
    if chain is None:
        raise UnsupportedStructure(
            "a branched substituent chain hanging off the ring is not "
            "supported yet"
        )

    chain_set = set(chain)
    chain_selenols = {s for s in selenols if next(iter(graph[s])) in chain_set}
    ring_selenols = selenols - chain_selenols
    if len(ring_selenols) < len(chain_selenols):
        ring_name, ring_is_compound = name_branch(
            graph, ring_atom, chain_root, {**halogens, **{s: "selanyl" for s in ring_selenols}}, mol=mol
        )
        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            se_locants = _se_locants(position_of, chain_selenols, graph)
            substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
            key, name = _candidate_key(chain_length, se_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    chain_name, chain_is_compound = name_branch(
        graph, chain_root, ring_atom, {**halogens, **{s: "selanyl" for s in chain_selenols}}, mol=mol
    )

    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            se_locants = _se_locants(position_of, ring_selenols, graph)
            substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
            key = _ring_candidate_key(ring_size, se_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_selenol(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            ring_selenols = _validate_and_collect_selenols(mol, aromatic_ring_atoms=ring_atoms)
            if len(ring_selenols) == 1:
                (only_se,) = ring_selenols
                (only_se_carbon,) = adjacency(mol)[only_se]
                if only_se_carbon in ring_atoms:
                    return _name_benzeneselenol(mol, ring_atoms)
            return _name_phenyl_chain_selenol(mol, ring_atoms)
    selenols = _validate_and_collect_selenols(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, _YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_eneselenol_carbons(graph, selenols, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        raise UnsupportedStructure(
            "polycyclic/spiro selenols are not supported yet (this module "
            "only handles acyclic chains and a single saturated ring)"
        )
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic selenol "
                "is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == _YNE_BOND_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a selenol is "
                "not supported yet -- only a ring double bond is in scope "
                "for this first pass (see P-31.1.3)"
            )
        ring_selenols = {s for s in selenols if next(iter(graph[s])) in ring_atoms}
        if not bonds and not ring_selenols:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than "
                    "the ring itself is not supported yet (see P-92)"
                )
            return _name_ring_substituent_chain_selenol(mol, selenols)
        if not bonds and ring_selenols and ring_selenols != selenols:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter alongside a ring-vs-chain selenol "
                    "comparison is not supported yet (see P-92)"
                )
            return _name_ring_with_selenol_chain_selenol(mol, selenols)
        if ring_selenols != selenols:
            raise UnsupportedStructure(
                "a selenol on a substituent branch chain rather than the "
                "ring itself is not supported yet"
            )
        return _name_cyclic_selenol(mol, selenols, stereo, bonds)

    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _se_locants(position_of, selenols, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _se_locants({a: i + 1 for i, a in enumerate(c)}, selenols, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every selenol-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain, "
            "or an -SeH expressed as a 'selanyl' substituent prefix, is "
            "not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            se_locants = _se_locants(position_of, selenols, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, selenols, mol=mol)
            key, name = _candidate_key(chain_length, se_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
