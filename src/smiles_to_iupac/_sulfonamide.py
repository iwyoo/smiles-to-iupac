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
  `_sulfonic_acid.py` apply (mononuclear parent, or a saturated two-carbon
  chain, regardless of other substituents), e.g. 'ethanesulfonamide'
  (PubChem structure match).
- P-44.4.1.8 / P-45.2: the -SO2NH2 locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as `_sulfonic_acid.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -SO2NH2 suffix. Confirmed via PubChem: 'propane-1-sulfonamide'
  (CCCS(=O)(=O)N), 'cyclohexanesulfonamide' (O=S(=O)(N)C1CCCCC1),
  '2-chloroethanesulfonamide' (ClCCS(=O)(=O)N).
- P-91.3/P-92: a molecule with
  one or more *specified* tetrahedral stereocenters -- every one on the
  principal chain/ring itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
  ascending locant order, same pattern as `_sulfonic_acid.py` (chain and
  ring both). The sulfonamide's own sulfur is never itself a potential
  stereocenter (its two =O substituents are identical), confirmed via
  RDKit `FindPotentialStereo` on `CC(C)S(=O)(=O)N`.
- The sulfonamide nitrogen may carry zero, one, or two plain,
  unsubstituted, saturated, acyclic alkyl substituents (branched or
  unbranched), each cited as its own 'N-'-prefixed substituent directly
  ahead of the parent stem, in alphanumerical order (P-14.5.2, ignoring
  italicized prefixes like 'tert-' via `alpha_sort_key`), with a 'di'
  multiplying prefix (and a single shared 'N,N-' pair) when both are
  identical -- exactly `_amide.py`'s own N-/N,N-disubstitution pattern
  (PR #335), built with `name_branch` (P-29 PIN style, fixed
  project-wide by PR #237). A *compound* N-substituent (has its own
  locant) is always parenthesized -- 'N-(propan-2-yl)methanesulfonamide',
  not PubChem's own raw 'N-propan-2-ylmethanesulfonamide' (CID 312702),
  same correction as `_urea.py` (see that module's docstring for the
  Blue Book citations). Confirmed via PubChem for the non-compound
  cases: 'N-methylmethanesulfonamide' (CS(=O)(=O)NC), 'N,N-
  dimethylmethanesulfonamide' (CS(=O)(=O)N(C)C),
  'N-tert-butylmethanesulfonamide' (CID 4130162),
  'N,N-di(propan-2-yl)methanesulfonamide' (CID 284325, a multiplied
  compound name parenthesized to avoid ambiguity),
  'N,N-ditert-butylmethanesulfonamide' (CID 58624041, a multiplied
  retained name not parenthesized),
  'N-tert-butyl-N-ethylmethanesulfonamide' (CID 58540473, alphabetized
  ignoring 'tert-').

Scope, deliberately narrow, mirroring `_sulfonic_acid.py`'s own first pass
exactly: a single -SO2NH2 on an acyclic chain or on a single saturated
carbon ring (monocyclic), with no other heteroatom anywhere in the molecule
except the sulfonamide group's own oxygens/nitrogen (and any N-alkyl
substituent's carbons) -- acid/amide-vs-other Table 3.3 seniority
coexistence is future work.

P-31.1.3: a monocyclic ring bearing a sulfonamide and exactly one C=C
ring double bond -- e.g. 'cyclohex-2-ene-1-sulfonamide',
'cyclohex-3-ene-1-sulfonamide', both confirmed via PubChem PUG REST.
Mirrors `_sulfonic_acid.py`'s identical extension; deliberately narrow: a
ring triple bond, any other ring substituent alongside the ring double
bond, and an N-substituent alongside a ring double bond (an unverified
combination) are all still explicitly rejected pending further
verification.

Explicitly out of scope (raise
`UnsupportedStructure`): an unsaturated or ring-bearing N-substituent (a
branched but otherwise plain saturated acyclic N-substituent is
supported, see above); an N-substituted sulfonamide, more than one
sulfonamide, or a specified stereocenter, on a von Baeyer polycyclic or
spiro skeleton (a single primary -SO2NH2 on such a skeleton is
supported, P-23/P-24 numbering integration via `_polycyclic_suffix.py`
with `elide_e=False`, mirroring `_thiol.py`); unsaturation reaching
outside the ring or a ring triple bond, a -SO2NH2 on a
substituent branch off an otherwise-unsubstituted *saturated* ring, two
or more -SO2NH2 groups, and a sulfonamide on a carbon that is also part
of a C=C/C#C bond. Two *aromatic*-ring cases:
`_name_phenyl_chain_sulfonamide` names a -SO2NH2 chain hanging off a
single plain, unsubstituted benzene ring (e.g.
'3-phenylpropane-1-sulfonamide'), mirroring `_sulfonic_acid.py`'s
identical benzene-ring-substituent path -- narrower than the acyclic
path: no N-alkyl substitution, no chain unsaturation, and no specified
stereocenter. `_name_benzenesulfonamide` names -SO2NH2 directly on a
benzene ring carbon (with or without other ring substituents), e.g.
'benzenesulfonamide' (PubChem PUG REST match for c1ccccc1S(=O)(=O)N),
'4-methylbenzenesulfonamide' (PUG REST match for Cc1ccc(cc1)S(=O)(=O)N),
mirroring `_sulfonic_acid.py`'s `_name_benzenesulfonic_acid` with the
retained name 'benzene' as stem -- the -SO2NH2's own locant is never
cited here; supports the same N-alkyl substitution as the acyclic path
(reusing `_n_alkyl_info`, merged into the ring citation via
`_add_n_names`), e.g.
'N-methylbenzenesulfonamide'/'N,N-dimethylbenzenesulfonamide' (both
PubChem PUG REST IUPACName matches) -- still no specified stereocenter.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    all_chains,
    bfs,
    bond_locant,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    lowest_locant_set,
    most_multiple_bonds,
    name_from_substituents,
    non_single_bonds,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    ring_chain_attachments_with_halogens,
    ring_hosting_anchors,
    separate_aromatic_monocycles,
    ring_cycle,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    suffix_body,
)
from ._bicyclic import find_bicyclic_core
from ._numerals import alkane_name, alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._substituents import (
    substituents_for_ring,
    branch_atom_locant,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
    ring_branch_stereo_display,
    substituents_for_chain,
)

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


def _validate_and_collect_sulfonamides(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_sulfonamide`'s benzene-
    ring-substituent path (see `_name_phenyl_chain_sulfonamide`) can
    reuse this same validation for the rest of the molecule. Empty by
    default, so every other caller's behavior is unchanged."""
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
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
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
    return sulfur.GetIdx(), carbon.GetIdx(), nitrogen.GetIdx(), n_alkyl_carbons


def _reject_enesulfonamide_carbon(graph, so2nh2_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if so2nh2_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a sulfonamide on a carbon that is also part of a C=C/C#C bond "
            "is out of scope for this module"
        )


def _n_alkyl_info(full_graph, full_carbon_graph, nitrogen_idx, n_alkyl_carbons, non_single_bond_atoms, halogens, mol=None):
    """Mirrors `_amide.py`'s own N-alkyl handling (PR #335): each
    N-substituent must be a plain, unsubstituted, saturated, acyclic
    alkyl chain (branched or unbranched). Returns (n_names,
    n_substituent_atoms) -- the latter must be excluded from the carbon
    graph before picking the principal chain, since these carbons hang
    off the (excluded) sulfonamide nitrogen rather than off the
    sulfonamide carbon itself."""
    n_names = []
    n_substituent_atoms = set()
    for n_alkyl_c in n_alkyl_carbons:
        n_atoms, _ = bfs(full_carbon_graph, n_alkyl_c)
        n_atoms = set(n_atoms)
        if n_atoms & non_single_bond_atoms:
            raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
        if any(nbr in n_atoms for h in halogens for nbr in full_graph[h]):
            # A halogen on the N-substituent is invisible to
            # `full_carbon_graph` -- calling `name_branch` with an empty
            # halogens dict would otherwise walk straight through it as
            # if it were a chain-extending atom, mirroring the bug
            # `_amide.py` found and fixed (PR #335).
            raise UnsupportedStructure(
                "a substituted N-substituent (e.g. bearing a halogen) is "
                "not supported yet; only a plain, unsubstituted alkyl "
                "N-substituent is in scope"
            )
        n_names.append(name_branch(full_graph, n_alkyl_c, nitrogen_idx, {}, mol=mol))
        n_substituent_atoms |= n_atoms
    return n_names, n_substituent_atoms


def _name_from_substituents(chain_length, so2nh2_locant, ene_locants, yne_locants, grouped, n_names=()):
    # Only chain_length == 1 omits a substituent prefix's own locant too
    # (see `_alcohol.py`'s equivalent comment) -- an 'N-' locant is never
    # omitted either way (see `_add_n_names`/`format_substituent_prefixes`:
    # it marks a different atom than the chain itself).
    prefix = format_substituent_prefixes(_add_n_names(grouped, n_names), omit_locants=chain_length == 1)
    return prefix + name_from_substituents(chain_length, ene_locants, yne_locants, "sulfonamide", [so2nh2_locant])


def _candidate_key(chain_length, so2nh2_locant, ene_locants, yne_locants, substituents, n_names=()):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, so2nh2_locant, ene_locants, yne_locants, grouped, n_names)
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

def _ring_name_from_substituents(ring_size, so2nh2_locant, ene_locants, yne_locants, grouped, n_names=()):
    has_unsaturation = bool(ene_locants or yne_locants)
    stem = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if not has_unsaturation:
        if total_subs == 0 and not n_names:
            # P-14.3.3: the sole substituent on an otherwise unsubstituted
            # ring has no locant to distinguish, e.g.
            # 'cyclohexanesulfonamide'.
            return stem + "sulfonamide"
        if total_subs == 0:
            prefix = format_substituent_prefixes(_add_n_names(grouped, n_names), omit_locants=True)
            return f"{prefix}{stem}sulfonamide"
        prefix = format_substituent_prefixes(_add_n_names(grouped, n_names))
        return f"{prefix}{stem}-{so2nh2_locant}-sulfonamide"

    # A competing ring double/triple bond (P-31.1.3) means the
    # sulfonamide's locant is never omittable even when it's the sole
    # substituent -- mirrors `_sulfonic_acid.py`'s identical treatment.
    unsaturated_stem = stem[:-3]
    needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    prefix = format_substituent_prefixes(_add_n_names(grouped, n_names))
    body = suffix_body(ene_locants, yne_locants, "sulfonamide", [so2nh2_locant])[0]
    return prefix + unsaturated_stem + ("a" if needs_stem_a else "") + "-" + body


def _ring_candidate_key(ring_size, so2nh2_locant, ene_locants, yne_locants, substituents, n_names=()):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, so2nh2_locant, ene_locants, yne_locants, grouped, n_names)
    return so2nh2_locant, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=None):
    return ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=mol)


def _name_cyclic_sulfonamide(mol, sulfur_idx, so2nh2_carbon, n_names, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid`),
    and the winning ring numbering's own locants for those atoms are used
    to format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_sulfonic_acid.py`'s own
    case): exactly one stereocenter on the ring's sole substituent branch
    instead embeds a bracketed descriptor into that substituent's own
    name, in place of the usual ring-locant prefix."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(substituents_for_ring(graph, ring_order, halogens, excluded, mol=mol).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "sulfonamide is not supported yet (see module docstring)"
        )
    if bonds and n_names:
        raise UnsupportedStructure(
            "an N-substituent alongside a ring double bond is not "
            "supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2nh2_locant = position_of[so2nh2_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, so2nh2_locant, ene_locants, yne_locants, substituents, n_names)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _add_n_names(grouped, n_names):
    """Merge N-alkyl substituent names into a *display-only* copy of the
    ring `grouped` dict, each getting the non-numeric locant 'N' --
    mirrors `_amine.py`'s identical `_add_n_names` (P-66.4's N-locant
    interleaving, PR #443). Ring-numbering ranking (`locant_set`/
    `citation_locants` below) stays on the original, N-less `grouped`,
    since an 'N' locant never affects which ring rotation wins -- only
    the returned name string needs it, e.g. 'N,4-dimethylbenzenesulfonamide'
    (PubChem PUG REST IUPACName match for `Cc1ccc(cc1)S(=O)(=O)NC`)."""
    if not n_names:
        return grouped
    display = {name: {"locants": list(info["locants"]), "compound": info["compound"]} for name, info in grouped.items()}
    for name, is_compound in n_names:
        info = display.setdefault(name, {"locants": [], "compound": is_compound})
        info["locants"].append("N")
    return display


def _benzenesulfonamide_name_from_substituents(grouped, n_names=()):
    # Mirrors `_sulfonic_acid.py`'s `_benzenesulfonic_acid_name_from_substituents`:
    # the mancude ring's own numbering is always free to start at the
    # -SO2NH2 carbon, so its locant is never cited even when other
    # substituents need theirs, e.g. '4-methylbenzenesulfonamide'.
    display = _add_n_names(grouped, n_names)
    if not display:
        return "benzenesulfonamide"
    return f"{format_substituent_prefixes(display)}benzenesulfonamide"


def _benzenesulfonamide_candidate_key(so2nh2_locant, substituents, n_names=()):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzenesulfonamide_name_from_substituents(grouped, n_names)
    return so2nh2_locant, locant_set, citation_locants, name


def _name_benzenesulfonamide(mol, ring_atoms, exempt_atoms=None):
    """P-65.3.1: -SO2NH2 attached directly to a benzene ring carbon -- e.g.
    'benzenesulfonamide', '4-methylbenzenesulfonamide', and now
    'N-methylbenzenesulfonamide'/'N,N-dimethylbenzenesulfonamide'/
    'N,4-dimethylbenzenesulfonamide' (all PubChem PUG REST IUPACName
    matches). Mirrors `_sulfonic_acid.py`'s `_name_benzenesulfonic_acid`
    with the retained name 'benzene' as stem; the -SO2NH2's own locant is
    never cited here. The N-alkyl substituent itself reuses
    `_n_alkyl_info` unchanged from the acyclic path (`name_sulfonamide`),
    but is merged into the same citation as any ring substituent (via
    `_add_n_names`) rather than a separate 'N-'/ring-prefix block, so an
    identically-named pair (e.g. two 'methyl's) collapses into one
    multiplied citation exactly like PubChem's own name."""
    sulfur_idx, so2nh2_carbon, nitrogen_idx, n_alkyl_carbons = _validate_and_collect_sulfonamides(
        mol, aromatic_ring_atoms=exempt_atoms or ring_atoms
    )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside benzenesulfonamide is not "
            "supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    non_single_bond_atoms = {a for a, b, _ in non_single_bonds(mol)} | {b for a, b, _ in non_single_bonds(mol)}
    n_names, n_substituent_atoms = _n_alkyl_info(
        graph, carbon_adjacency(mol), nitrogen_idx, n_alkyl_carbons, non_single_bond_atoms, halogens, mol=mol
    )
    excluded = {sulfur_idx} | n_substituent_atoms
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2nh2_locant = position_of[so2nh2_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _benzenesulfonamide_candidate_key(so2nh2_locant, substituents, n_names)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_sulfonamide(mol, ring_atoms):
    """Name a sulfonamide whose -SO2NH2 lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropane-1-sulfonamide.
    The ring is cited as a 'phenyl' substituent prefix on the chain,
    which is the parent hydride, mirroring
    `_sulfonic_acid.py`'s `_name_phenyl_chain_sulfonic_acid`. Narrower
    than the acyclic path above: no N-alkyl substitution, no chain
    unsaturation, and no specified stereocenter -- each is a separate
    follow-up (see
    tasks/phenyl-substituent-on-sulfonamide-chain.md's scope note)."""
    sulfur_idx, so2nh2_carbon, _, n_alkyl_carbons = _validate_and_collect_sulfonamides(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if n_alkyl_carbons:
        raise UnsupportedStructure(
            "an N-alkyl-substituted sulfonamide alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "sulfonamide chain is not supported yet"
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
            "sulfonamide chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_chain_attachments_with_halogens(graph, rings, set(), halogens)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain sulfonamide is not "
            "supported yet"
        )
    chain, branches = longest_branched_chain_through(graph, so2nh2_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        so2nh2_locant = position_of[so2nh2_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, so2nh2_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_sulfonamide(mol, sulfur_idx, so2nh2_carbon):
    """Name a sulfonamide whose -SO2NH2 lies entirely on a single branched
    chain hanging off one atom of an otherwise-plain saturated monocyclic
    ring (the ring itself bears no sulfonamide) -- e.g.
    cyclohexylmethanesulfonamide. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_sulfonamide` above and
    `_sulfonic_acid.py`'s `_name_ring_substituent_chain_sulfonic_acid`.
    Narrower than the benzene-ring case: no N-alkyl substitution."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])
    excluded = {sulfur_idx}

    attachment = ring_chain_attachment(graph, ring_atoms, excluded)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment

    chain, branches = longest_branched_chain_through(graph, so2nh2_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
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
        so2nh2_locant = position_of[so2nh2_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, so2nh2_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_von_baeyer_or_spiro_sulfonamide(mol, sulfur_idx, so2nh2_carbon, n_alkyl_carbons, bonds, stereo):
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single sulfonamide (-SO2NH2) suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`'s
    `elide_e=False` ('sulfonamide' begins with a consonant, P-16.3.3,
    same as `_thiol.py`'s 'thiol' -- e.g. 'bicyclo[2.2.1]heptane-2-
    sulfonamide', PubChem CID 45080580). Mirrors `_sulfinamide.py`'s own
    bicyclic/polycyclic-before-spiro dispatch order and restrictions:
    exactly one primary sulfonamide on the ring system itself; ring
    unsaturation composes on the bicyclic/polycyclic branch (mirrors
    `_ketone.py`'s identical extension) but not the monospiro branch,
    which has no base mechanism yet (P-31.1.5)."""
    if n_alkyl_carbons:
        raise UnsupportedStructure(
            "an N-substituted sulfonamide on a von Baeyer bicyclic/"
            "polycyclic or monospiro ring system is not supported yet"
        )

    bicyclic_core = find_bicyclic_core(mol)
    polycyclic_core = None
    von_baeyer_ring_count = None
    if bicyclic_core is None:
        for candidate_ring_count in (3, 4, 5, 6):
            polycyclic_core = find_polycyclic_core(mol, candidate_ring_count)
            if polycyclic_core is not None:
                von_baeyer_ring_count = candidate_ring_count
                break
    if bicyclic_core is not None or polycyclic_core is not None:
        return name_von_baeyer_suffix(
            mol,
            so2nh2_carbon,
            {sulfur_idx},
            "sulfonamide",
            "sulfonamide",
            bicyclic_core,
            polycyclic_core,
            von_baeyer_ring_count,
            elide_e=False,
            stereo=stereo,
            bonds=bonds,
        )

    if bonds:
        raise UnsupportedStructure(
            "an unsaturated monospiro ring system is not supported yet "
            "(see P-31.1.5)"
        )

    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return name_monospiro_suffix(
            mol, so2nh2_carbon, {sulfur_idx}, "sulfonamide", "sulfonamide", spiro_atom, elide_e=False, stereo=stereo
        )

    raise UnsupportedStructure(
        "polycyclic and fused-ring sulfonamides are not supported yet "
        "(P-23/P-25 numbering integration with a suffix group is future "
        "work)"
    )


def name_sulfonamide(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        union = set().union(*aromatic_rings)
        sulfur_idx, *_ = _validate_and_collect_sulfonamides(mol, aromatic_ring_atoms=union)
        anchors = [sulfur_idx]
        host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
        if host is not None:
            return _name_benzenesulfonamide(mol, host, union)
        return _name_phenyl_chain_sulfonamide(mol, union)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            _, so2nh2_carbon, _, _ = _validate_and_collect_sulfonamides(mol, aromatic_ring_atoms=ring_atoms)
            if so2nh2_carbon in ring_atoms:
                return _name_benzenesulfonamide(mol, ring_atoms)
            return _name_phenyl_chain_sulfonamide(mol, ring_atoms)
    sulfur_idx, so2nh2_carbon, nitrogen_idx, n_alkyl_carbons = _validate_and_collect_sulfonamides(mol)
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
    halogens = halogen_substituents(mol)
    n_names, n_substituent_atoms = _n_alkyl_info(
        graph, carbon_adjacency(mol), nitrogen_idx, n_alkyl_carbons, non_single_bond_atoms, halogens, mol=mol
    )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        return _name_von_baeyer_or_spiro_sulfonamide(
            mol, sulfur_idx, so2nh2_carbon, n_alkyl_carbons, bonds, stereo
        )
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
        if any(order == _YNE_ORDER for _, _, order in ring_bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a sulfonamide "
                "is not supported yet -- only a ring double bond is in "
                "scope for this first pass (see P-31.1.3)"
            )
        if so2nh2_carbon not in ring_atoms:
            if not bonds and not n_alkyl_carbons:
                if stereo is not None:
                    raise UnsupportedStructure(
                        "a stereocenter on a substituent branch rather "
                        "than the ring itself is not supported yet (see "
                        "P-92)"
                    )
                return _name_ring_substituent_chain_sulfonamide(mol, sulfur_idx, so2nh2_carbon)
            raise UnsupportedStructure(
                "a sulfonamide on a substituent branch chain rather than "
                "the ring itself is not supported yet"
            )
        return _name_cyclic_sulfonamide(mol, sulfur_idx, so2nh2_carbon, n_names, stereo, ring_bonds)

    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    full_carbon_graph = carbon_adjacency(mol)
    carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in n_substituent_atoms}
    chains = all_chains(carbon_graph)
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if so2nh2_carbon not in chain:
            continue
        if stereo is not None and any(atom not in chain for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            so2nh2_carbon in c for c in chains
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
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2nh2_locant = position_of[so2nh2_carbon]
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, so2nh2_locant, ene_locants, yne_locants, substituents, n_names)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
