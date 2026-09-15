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
amine whose extra N-substituent(s) are saturated alkyl chains, branched or
straight, halogen-substituted or not: the N-linked carbon starting the
largest carbon skeleton becomes the parent chain (suffixed '-amine' as
usual), and each other N-linked chain is named the same way an ordinary
substituent branch is (`name_branch`) and cited as an 'N-'/'N,N-'
substituent prefix, e.g. 'N-ethylethanamine' (diethylamine,
PubChem-verified), 'N,N-dimethylmethanamine' (trimethylamine,
PubChem-verified), 'N-(propan-2-yl)propan-1-amine' (PubChem CID 89119, raw
'N-propan-2-ylpropan-1-amine' corrected to this project's established
compound-prefix parenthesization, see `_amide.py`), and
'N-(2-chloroethyl)propan-1-amine' (PubChem CID 3045065). An N-substituent
may also carry a single C=C/C#C bond as long as it's otherwise unbranched
(`name_branch` has no ene/yne machinery of its own yet, so this narrow
case is named directly by `unbranched_unsaturated_substituent_name`
instead), e.g. 'N-prop-2-enylpentan-1-amine' (PubChem CID 12442641),
'N-prop-2-ynylbutan-1-amine' (PubChem CID 3465926). A halogen substituent
on the parent chain also coexists with this (P-14.5.2: the 'N-' prefix
interleaves alphabetically with any other substituent prefix rather than
always citing first, e.g. 'ClCCNCC' -> '2-chloro-N-ethylethanamine',
PubChem-verified) -- `_add_n_names`/`format_substituent_prefixes`'s
non-numeric 'N' locant handles this generically, so it applies equally to
`_ammonium.py`'s reuse of this same machinery. It otherwise mirrors
`_alcohol.py`'s scope restrictions:

Explicitly out of scope (raise `UnsupportedStructure`):
- A secondary/tertiary amine nitrogen with a *branched* unsaturated
  N-substituent, more than one multiple bond on the same N-substituent,
  an N-substituent nitrogen directly on the alkene/alkyne carbon itself
  (an enamine-shaped nitrogen, rejected by `_reject_enamine_carbons`
  regardless of which side of the amine it's on), or more than two
  N-substituents (mirrors `_amide.py`'s own N-substituent restrictions).
- A secondary/tertiary amine nitrogen on or attached to a ring, or
  coexisting with another amine nitrogen elsewhere in the molecule (a
  diamine where one nitrogen is secondary/tertiary) — both deferred as
  separate, larger extensions.
- Any other heteroatom (O, S, ...) — including molecules that would also
  need a senior characteristic group (Table 3.3); this module rejects those
  outright rather than attempting suffix-vs-suffix seniority.
- -NHR/-NR2 (N-substituted) on an aromatic ring, or a second amine
  nitrogen alongside a ring one — a plain, unsubstituted -NH2 directly on
  an otherwise-plain benzene ring is handled by `_name_aniline`
  (P-62.2.1.1.1), a separate follow-up from this narrower N-substituted/
  multi-group combination.
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

P-31.1.3: a monocyclic ring bearing a *primary* amine and exactly one C=C
ring double bond -- e.g. 'cyclohex-2-en-1-amine', 'cyclohex-3-en-1-amine',
both confirmed via PubChem PUG REST. Mirrors `_ketone.py`'s identical
extension; this is a narrow, separate slice from the existing
secondary/tertiary-amine-on-a-ring restriction above (which still applies
unchanged -- this axis only opens up for the plain -NH2 case).
Deliberately narrow: a ring triple bond, and any other substituent
alongside the ring double bond, are both still explicitly rejected
pending further verification.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locant,
    bond_locants,
    carbon_adjacency,
    component_subgraph,
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
    ring_chain_attachment_with_halogens,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    suffix_body,
    unbranched_unsaturated_substituent_name,
)
from ._numerals import alkane_name, alkyl_name
from ._substituents import (
    branch_atom_locant,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
    ring_branch_stereo_display,
)

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _validate_and_collect_amines(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (amines, n_carbons_by_nitrogen): the set of amine-nitrogen
    atom indices, and a dict mapping each to a tuple of its 1-3 carbon
    neighbor indices (the -NH2 carbon for a primary amine; the parent-chain
    carbon plus 0-2 N-substituent carbons for a secondary/tertiary one).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_amine`'s benzene-ring-substituent
    path (see `_name_phenyl_chain_amine`) can reuse this same validation for
    the rest of the molecule. Empty by default, so every other caller's
    behavior is unchanged. Mirrors `_thiol.py`'s
    `_validate_and_collect_thiols`."""
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
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
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


def _add_n_names(grouped, n_names, n_locants=None):
    """Merge N-substituent names (P-66.4) into a parent-chain `grouped`
    dict, each getting a non-numeric locant -- for display only
    (`format_substituent_prefixes`'s alphabetical citation, P-14.5.2, sorts
    by name regardless of locant type, and `_locant_sort_key` in
    `_substituents.py` puts a non-numeric locant after every numeric one
    when the same name also occurs on the numbered chain). Never mutates
    the input `grouped` -- callers still need the original, numeric-only
    version for both chain-orientation ranking and the `chain_length`-based
    locant-omission branches below, neither of which this locant may enter
    (it doesn't move when the chain is renumbered, and its presence
    doesn't make the amine's own carbon locant any less omittable).

    `n_names`: (name, is_compound) pairs, one per N-substituent -- the same
    shape `name_branch` returns, so a branched/compound N-substituent (e.g.
    'propan-2-yl') gets its own enclosing marks via
    `format_substituent_prefixes` exactly like any other compound prefix.

    `n_locants`: parallel to `n_names`, or None. A single amine nitrogen
    (the only kind this project supported until P-16.9.2's two-amine case)
    has no locant of its own to disambiguate, so its N-substituents just
    get the bare 'N' -- `n_locants=None`, or an entry of `None`, both mean
    that. With two coexisting amine nitrogens, each entry is instead that
    nitrogen's own parent-chain locant, producing 'N<k>' (e.g. 'N2') per
    P-16.9.2 (superseding the older N/N'-prime convention for this exact
    case -- see `_name_multi_amine_chain`'s docstring)."""
    if not n_names:
        return grouped
    display = {name: {"locants": list(info["locants"]), "compound": info["compound"]} for name, info in grouped.items()}
    for i, (name, is_compound) in enumerate(n_names):
        locant = n_locants[i] if n_locants is not None else None
        display_locant = "N" if locant is None else f"N{locant}"
        info = display.setdefault(name, {"locants": [], "compound": is_compound})
        info["locants"].append(display_locant)
    return display


def _name_from_substituents(chain_length, amine_locants, ene_locants, yne_locants, grouped, n_names=(), n_locants=None):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants (prefix or suffix)
        # are always '1' and never cited -- an 'N-' locant is never
        # omittable though (see `_add_n_names`/`format_substituent_prefixes`,
        # it marks a different atom than the mononuclear carbon itself).
        amine_word = multiplied_word(len(amine_locants), "amine")
        stem = alkane_name(1)
        if amine_word[0] in "aeiouy":
            stem = stem[:-1]
        prefix = format_substituent_prefixes(_add_n_names(grouped, n_names, n_locants), omit_locants=True)
        return prefix + stem + amine_word

    if chain_length == 2 and not has_unsaturation and total_subs == 0 and len(amine_locants) == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly one
        # substituent (here, the sole -NH2) in total has only one possible
        # structure, so the locant is omittable, e.g. 'ethanamine' -- an
        # N-substituent doesn't change this (it doesn't create any
        # carbon-position asymmetry), but still needs its own 'N-' prefix,
        # e.g. 'N-ethylethanamine'.
        prefix = format_substituent_prefixes(_add_n_names(grouped, n_names, n_locants))
        return prefix + alkane_name(2)[:-1] + "amine"

    prefix = format_substituent_prefixes(_add_n_names(grouped, n_names, n_locants))
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = suffix_body(ene_locants, yne_locants, multiplied_word(len(amine_locants), "amine"), amine_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, amine_locants, ene_locants, yne_locants, substituents, n_names=(), n_locants=None):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first.

    `n_names`/`n_locants`: N-substituent names (P-66.4) from a secondary/
    tertiary amine's other N-linked chain(s) (`n_locants` parallel, or None
    for a single amine nitrogen -- see `_add_n_names`), passed straight
    through to `_name_from_substituents` (which folds them into the
    *displayed* name via `_add_n_names` so they interleave alphabetically
    with any parent-chain substituent prefix, P-14.5.2, e.g.
    '2-chloro-N-ethyl...') but excluded from every locant-set computation
    below since an N-locant is never a candidate for the lowest-locant
    chain-orientation tie-break."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    amine_locant_set = lowest_locant_set(amine_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, amine_locants, ene_locants, yne_locants, grouped, n_names, n_locants)
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


def _substituents_for_chain(graph, chain, halogens, amines, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in amines]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _best_chain_name(
    carbon_graph, graph, halogens, amines, bonds, stereo=None, n_names=(), mol=None, n_names_by_nitrogen=None
):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch is out of scope). Returns (name, position_of) so a
    caller wrapping a stereo prefix around `name` can still place it
    outermost using the winning chain's own locants.

    `n_names`: a secondary/tertiary amine's other N-substituent name(s)
    (P-66.4), passed straight through to `_candidate_key` so they're
    already correctly interleaved into `name` itself (P-14.5.2) -- the
    caller no longer glues an 'N-' prefix on afterward. Ignored when
    `n_names_by_nitrogen` is given.

    `n_names_by_nitrogen`: for `_name_multi_amine_chain`'s coexisting-
    amines case (P-16.9.2) -- `{nitrogen_atom_idx: [(name, is_compound),
    ...]}`. Each nitrogen's own chain locant depends on which orientation
    of which candidate chain wins, so (unlike the single-amine `n_names`
    above) this can't be resolved by the caller ahead of time -- it's
    recomputed per candidate below, from that candidate's own
    `position_of`, into the flat `n_names`/`n_locants` pair
    `_candidate_key` expects."""
    chains = longest_chains(carbon_graph)
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _amine_locants(position_of, amines, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
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
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, amines, mol=mol)
            if n_names_by_nitrogen is not None:
                candidate_n_names = []
                candidate_n_locants = []
                for n_idx, locant in zip(amines, amine_locants):
                    for sub_name, is_compound in n_names_by_nitrogen.get(n_idx, []):
                        candidate_n_names.append((sub_name, is_compound))
                        candidate_n_locants.append(locant)
            else:
                candidate_n_names, candidate_n_locants = n_names, None
            key, name = _candidate_key(
                chain_length, amine_locants, ene_locants, yne_locants, substituents, candidate_n_names, candidate_n_locants
            )
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of
    return best_name, best_position_of


def _name_acyclic_secondary_tertiary_amine(mol, n_idx, n_carbons, bonds, stereo=None):
    """Name a secondary/tertiary amine: the N-linked carbon starting the
    largest carbon skeleton becomes the parent chain (suffixed '-amine' via
    `_best_chain_name`, same as a primary amine), and each other N-linked
    chain -- a branched and/or halogen-substituted saturated alkyl, named
    via `name_branch` the same way an ordinary substituent branch is
    (P-66.4/P-29) -- is cited as an 'N-'/'N,N-' prefix, e.g.
    'N-ethylethanamine' (diethylamine), 'N-(propan-2-yl)propan-1-amine'
    (PubChem CID 89119, raw 'N-propan-2-ylpropan-1-amine' corrected to this
    project's established compound-prefix parenthesization, see
    `_amide.py`), 'N-(2-chloroethyl)propan-1-amine' (PubChem CID 3045065).
    An unbranched N-substituent may also carry one C=C/C#C bond, named
    directly by `unbranched_unsaturated_substituent_name` instead
    (`name_branch` has no ene/yne machinery of its own yet) -- e.g.
    'N-prop-2-enylpentan-1-amine' (PubChem CID 12442641).

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
    components = {c: component_subgraph(full_carbon_graph, c) for c in n_carbons}
    parent_root = max(n_carbons, key=lambda c: len(components[c]))
    other_roots = [c for c in n_carbons if c != parent_root]

    n_names = []
    for other in other_roots:
        other_atoms = set(components[other])
        other_bonds = [b for b in non_single_bonds(mol) if b[0] in other_atoms and b[1] in other_atoms]
        if other_bonds:
            if len(other_bonds) > 1:
                raise UnsupportedStructure(
                    "more than one multiple bond on an N-substituent is "
                    "not supported yet"
                )
            name = unbranched_unsaturated_substituent_name(full_carbon_graph, other, other_bonds[0])
            if name is None:
                raise UnsupportedStructure(
                    "a branched unsaturated N-substituent is not supported yet"
                )
            n_names.append((name, False))
            continue
        n_names.append(name_branch(graph, other, n_idx, halogens, mol=mol))

    excluded_atoms = set()
    for other in other_roots:
        excluded_atoms |= set(components[other])
    parent_carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in excluded_atoms}
    # A multiple bond entirely on an N-substituent's own component (now
    # named directly by `unbranched_unsaturated_substituent_name`
    # above) isn't part of the parent chain at all -- passing it through
    # unfiltered would make `_best_chain_name`'s `bond_locants` unable to
    # place it on any candidate chain, since none of them contain its
    # atoms, incorrectly rejecting every chain as ineligible.
    parent_bonds = [b for b in bonds if b[0] in parent_carbon_graph and b[1] in parent_carbon_graph]

    best_name, best_position_of = _best_chain_name(
        parent_carbon_graph, graph, halogens, {n_idx}, parent_bonds, stereo, n_names=n_names, mol=mol
    )

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_multi_amine_chain(mol, amines, n_carbons_by_nitrogen, bonds, stereo=None):
    """Two or more coexisting amine nitrogens on one shared parent chain,
    at least one of which is secondary/tertiary (e.g. 'NCCN(C)C',
    N2,N2-dimethylethane-1,2-diamine; 'NCC(N)C(N)CN(C)C',
    N1,N1-dimethylbutane-1,2,3,4-tetramine). Each extra N-substituent is
    named exactly like `_name_acyclic_secondary_tertiary_amine` already
    does for a single amine -- via `name_branch` on its own connected-
    component root, isolated from the rest of the carbon skeleton by the
    same "N isn't in the carbon-only graph" property that function relies
    on -- but is now tagged with *that nitrogen's own chain locant* rather
    than a bare 'N', per P-16.9.2:

    > "Superscript arabic numbers, which are the locants of the parent
    > structure, are used to differentiate the nitrogen atoms of di- and
    > polyamines... except for geminal amines." Worked PIN example:
    > 'N1-ethyl-N2-methylethane-1,2-diamine' -- superseding the older
    > N/N'-prime convention (still seen in some database autonames, e.g.
    > PubChem's 'N,N'-dimethyl...' for the same shape) for this exact
    > case.

    A nitrogen's locant depends on which candidate chain orientation wins
    (P-44 lowest-locant tie-breaks), so it can't be resolved here --
    `_best_chain_name`'s `n_names_by_nitrogen` recomputes it per candidate
    instead of this function gluing on a fixed prefix like the
    single-amine case does.

    Scope, deliberately narrow (see `tasks/amine-cross-cutting-
    generalization.md`): every amine nitrogen must contribute to one
    shared chain (geminal -- two nitrogens on the same carbon -- is
    explicitly excluded by P-16.9.2 itself and raises here, as does a
    nitrogen not reachable from every other one via a single carbon
    backbone); each nitrogen's extra substituents are plain saturated
    (branched/halogenated) alkyl only, no unsaturation (mirroring
    `_name_acyclic_secondary_tertiary_amine`'s own unbranched-only
    unsaturated-N-substituent support would need its own per-nitrogen
    locant plumbing there too -- future work)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    full_carbon_graph = carbon_adjacency(mol)

    all_n_carbons = [c for n in amines for c in n_carbons_by_nitrogen[n]]
    components = {c: component_subgraph(full_carbon_graph, c) for c in all_n_carbons}
    # The shared parent chain is whichever connected carbon component has a
    # neighbor from *every* nitrogen -- not simply the largest component,
    # which a same-sized N-substituent fragment (e.g. an N-ethyl tied with
    # a 2-carbon real backbone) could otherwise be mistaken for.
    distinct_components = {frozenset(comp) for comp in components.values()}
    shared_components = [
        comp
        for comp in distinct_components
        if all(any(c in comp for c in n_carbons_by_nitrogen[n]) for n in amines)
    ]
    if len(shared_components) != 1:
        raise UnsupportedStructure(
            "every amine nitrogen must share exactly one connected carbon "
            "backbone for this naming path"
        )
    chain_component = shared_components[0]

    n_names_by_nitrogen = {}
    excluded_atoms = set()
    chain_anchors = set()
    for n in amines:
        on_chain = [c for c in n_carbons_by_nitrogen[n] if c in chain_component]
        if len(on_chain) != 1:
            raise UnsupportedStructure(
                "two or more amine nitrogens sharing the same carbon (geminal) is "
                "excluded from P-16.9.2's superscript-locant convention "
                "and is not supported yet"
            )
        chain_anchors.add(on_chain[0])
        extra_roots = [c for c in n_carbons_by_nitrogen[n] if c not in chain_component]
        names = []
        for root in extra_roots:
            root_component = set(components[root])
            if any(a in root_component and b in root_component for a, b, _ in bonds):
                raise UnsupportedStructure(
                    "an unsaturated N-substituent alongside another coexisting amine "
                    "nitrogen is not supported yet"
                )
            names.append(name_branch(graph, root, n, halogens, mol=mol))
            excluded_atoms |= root_component
        n_names_by_nitrogen[n] = names
    if len(chain_anchors) != len(amines):
        raise UnsupportedStructure(
            "two or more amine nitrogens sharing the same carbon (geminal) is "
            "excluded from P-16.9.2's superscript-locant convention "
            "and is not supported yet"
        )

    parent_carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in excluded_atoms}
    parent_bonds = [b for b in bonds if b[0] in parent_carbon_graph and b[1] in parent_carbon_graph]

    best_name, best_position_of = _best_chain_name(
        parent_carbon_graph, graph, halogens, amines, parent_bonds, stereo,
        mol=mol, n_names_by_nitrogen=n_names_by_nitrogen,
    )
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

    if len(amines) >= 2 and any(len(n_carbons_by_nitrogen[n]) > 1 for n in amines):
        return _name_multi_amine_chain(mol, amines, n_carbons_by_nitrogen, bonds, stereo)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    best_name, best_position_of = _best_chain_name(carbon_adjacency(mol), graph, halogens, amines, bonds, stereo, mol=mol)
    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _substituents_for_ring(graph, ring_order, halogens, amines, mol=None):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in amines]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, amine_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size,
        ene_locants,
        yne_locants,
        prefix,
        total_subs,
        multiplied_word(len(amine_locants), "amine"),
        amine_locants,
    )


def _ring_candidate_key(ring_size, amine_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    amine_locant_set = lowest_locant_set(amine_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, amine_locants, ene_locants, yne_locants, grouped)
    return amine_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, amines, stereo, halogens, mol=None):
    return ring_branch_stereo_display(graph, ring_order, amines, stereo, halogens, mol=mol)


def _name_cyclic_amine(mol, amines, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid`),
    and the winning ring numbering's own locants for those atoms are used
    to format a "(<locant><R/S>,...)-" prefix onto the name. The one
    narrow exception (`_ring_branch_stereo_display`, mirroring
    `_alcohol.py`/`_ketone.py`/`_thiol.py`/`_sulfonic_acid.py`'s own case):
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
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, amines, stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(_substituents_for_ring(graph, ring_order, halogens, amines, mol=mol).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "primary amine is not supported yet (see module docstring)"
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
            substituents = _substituents_for_ring(graph, candidate, halogens, amines, mol=mol)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, amine_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _aniline_name_from_substituents(grouped):
    # P-62.2.1.1.1: the retained name 'aniline' stands for the whole
    # ring+NH2 system (like 'phenol'/'phenoxide' in `_alcohol.py`/
    # `_alkoxide.py`), so the amine's own ring locant is never cited --
    # e.g. '4-methylaniline (PIN)', not '4-methylaniline-1-amine'.
    if not grouped:
        return "aniline"
    return f"{format_substituent_prefixes(grouped)}aniline"


def _aniline_candidate_key(amine_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _aniline_name_from_substituents(grouped)
    return amine_locant, locant_set, citation_locants, name


def _name_aniline(mol, ring_atoms):
    """P-62.2.1.1.1: -NH2 attached directly to a benzene ring carbon --
    e.g. 'aniline' (PubChem CID 6115), '4-methylaniline' (CID 7864), '2-
    chloroaniline' (CID 8863). Mirrors `_alcohol.py`'s `_name_phenol`
    exactly, with the retained name 'aniline' standing in for 'phenol'.
    Narrower than the general ring case: exactly one *primary* amine
    (directly on the ring, no N-alkyl substituent -- N-substituted
    aniline, e.g. 'N-methylaniline (PIN)', is a separate follow-up), and
    no specified stereocenter."""
    amines, n_carbons_by_nitrogen = _validate_and_collect_amines(mol, aromatic_ring_atoms=ring_atoms)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "more than one amine nitrogen directly on a benzene ring is "
            "not supported yet"
        )
    (n_idx,) = amines
    if len(n_carbons_by_nitrogen[n_idx]) > 1:
        raise UnsupportedStructure(
            "an N-substituted aniline is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure("a specified stereocenter alongside aniline is not supported yet")

    (amine_carbon,) = n_carbons_by_nitrogen[n_idx]
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {n_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            amine_locant = position_of[amine_carbon]
            substituents = _substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _aniline_candidate_key(amine_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_ring_substituent_chain_amine(mol, amines, n_carbons_by_nitrogen):
    """Name a primary amine whose -NH2 group(s) lie entirely on a single
    branched chain hanging off one atom of an otherwise-plain saturated
    monocyclic ring (the ring itself bears no amine) -- e.g.
    cyclohexylmethanamine. The ring is cited as a "cyclo..." substituent
    prefix (P-29.3.3) on the chain, which is the parent hydride, mirroring
    `_alcohol.py`'s `_name_ring_substituent_chain_alcohol`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, amines)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    anchor_n = next(iter(amines))
    (anchor_carbon,) = n_carbons_by_nitrogen[anchor_n]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, amines, halogens=halogen_substituents(mol))
    chain_set = set(chain)
    for n in amines:
        (carbon,) = n_carbons_by_nitrogen[n]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "an amine outside the single branched chain hanging off "
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
        amine_locants = _amine_locants(position_of, amines, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, amine_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_with_amine_chain_amine(mol, amines, n_carbons_by_nitrogen):
    """Name a primary-amine compound where the ring itself bears at least
    as many amines as a single unbranched chain hanging off exactly one
    ring atom does (P-44.1.1: the candidate with the greater count of the
    principal characteristic group amine is senior; P-44.1.2.2 resolves
    an exact tie in the ring's favor) -- e.g.
    2-(aminomethyl)cyclohexan-1-amine, 1-(aminomethyl)cyclohexane-1,2-
    diamine. The ring is the parent; the chain is cited as an
    '(amino...alkyl)' substituent prefix, mirroring `_alcohol.py`'s
    `_name_ring_with_hydroxy_chain_alcohol`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, amines)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, amines)
    if chain is None:
        raise UnsupportedStructure(
            "a branched substituent chain hanging off the ring is not "
            "supported yet"
        )

    chain_set = set(chain)
    chain_amines = {n for n in amines if next(iter(n_carbons_by_nitrogen[n])) in chain_set}
    ring_amines = amines - chain_amines
    if len(ring_amines) < len(chain_amines):
        # P-44.1.1: the chain captures strictly more amines, so it's the
        # senior parent and the ring (with its own one or more amines) is
        # cited as a substituent instead -- mirrors
        # `_name_ring_substituent_chain_amine` exactly, substituting the
        # ring's own name_branch-computed name for the plain "cyclo..."
        # one that function uses.
        ring_name, ring_is_compound = name_branch(
            graph, ring_atom, chain_root, {**halogens, **{n: "amino" for n in ring_amines}}, mol=mol
        )
        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            amine_locants = _amine_locants(position_of, chain_amines, graph)
            substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
            key, name = _candidate_key(chain_length, amine_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    chain_name, chain_is_compound = name_branch(
        graph, chain_root, ring_atom, {**halogens, **{n: "amino" for n in chain_amines}}, mol=mol
    )

    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            amine_locants = _amine_locants(position_of, ring_amines, graph)
            substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
            key = _ring_candidate_key(ring_size, amine_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_amine(mol, ring_atoms):
    """Name a primary amine whose -NH2 lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropan-1-amine. The ring is
    cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_thiol.py`'s `_name_phenyl_chain_thiol`. Narrower than the
    acyclic path above: exactly one *primary* amine (no N-alkyl
    substituent -- a secondary/tertiary amine alongside a ring is out of
    scope, same as `name_amine`'s existing ring guard below), no chain
    unsaturation, and no specified stereocenter."""
    amines, n_carbons_by_nitrogen = _validate_and_collect_amines(mol, aromatic_ring_atoms=ring_atoms)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "more than one amine nitrogen alongside a benzene-ring "
            "substituent is not supported yet"
        )
    (n_idx,) = amines
    if len(n_carbons_by_nitrogen[n_idx]) > 1:
        raise UnsupportedStructure(
            "a secondary/tertiary amine nitrogen on or attached to a ring "
            "is out of scope for this module"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "amine chain is not supported yet"
        )
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "amine chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain amine is not "
            "supported yet"
        )
    (carbon,) = n_carbons_by_nitrogen[n_idx]
    chain, branches = longest_branched_chain_through(graph, carbon, ring_atoms, amines, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        amine_locants = _amine_locants(position_of, amines, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, amine_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_amine(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            amines, n_carbons_by_nitrogen = _validate_and_collect_amines(mol, aromatic_ring_atoms=ring_atoms)
            if len(amines) == 1:
                (n_idx,) = amines
                if len(n_carbons_by_nitrogen[n_idx]) == 1:
                    (only_carbon,) = n_carbons_by_nitrogen[n_idx]
                    if only_carbon in ring_atoms:
                        return _name_aniline(mol, ring_atoms)
            return _name_phenyl_chain_amine(mol, ring_atoms)
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
    if num_rings == 0:
        return _name_acyclic_amine(mol, amines, n_carbons_by_nitrogen, bonds, stereo)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic amine "
                "is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == _YNE_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside an amine is "
                "not supported yet -- only a ring double bond is in scope "
                "for this first pass (see P-31.1.3)"
            )
        ring_amines = {
            n for n in amines if next(iter(n_carbons_by_nitrogen[n])) in ring_atoms
        }
        if not bonds and not ring_amines:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than "
                    "the ring itself is not supported yet (see P-92)"
                )
            return _name_ring_substituent_chain_amine(mol, amines, n_carbons_by_nitrogen)
        if not bonds and ring_amines and ring_amines != amines:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter alongside a ring-vs-chain amine "
                    "comparison is not supported yet (see P-92)"
                )
            return _name_ring_with_amine_chain_amine(mol, amines, n_carbons_by_nitrogen)
        return _name_cyclic_amine(mol, amines, stereo, bonds)
    raise UnsupportedStructure(
        "polycyclic and spiro amines are not supported yet (P-23/P-24/P-25 "
        "numbering integration with a suffix group is future work)"
    )
