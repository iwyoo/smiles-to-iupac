"""Naming of sulfinamides (the '-sulfinamide' suffix, -S(=O)NH2) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.3.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf), the
  nitrogen analogue of the sulfinic acid entry in the same table:
  'sulfinamide' is the preselected suffix for -S(=O)-NH2, one oxidation
  state below sulfonamide (`_sulfonamide.py`) -- the sulfur carries one
  carbon, one double-bonded oxygen, and one primary amide nitrogen (degree
  3, not 4). This module doesn't implement acid/amide-vs-other Table 3.3
  seniority yet (see scope note below), mirroring `_sulfinic_acid.py`'s
  own deferral.
- Like 'sulfinic acid'/'sulfonamide', 'sulfinamide' is cited as the parent
  hydride name followed directly by the suffix with no elision -- 'methane'
  + 'sulfinamide' -> 'methanesulfinamide' (PubChem structure match).
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_sulfinic_acid.py`/`_sulfonamide.py` apply (mononuclear parent, or a
  saturated two-carbon chain, regardless of other substituents), e.g.
  'ethanesulfinamide' (PubChem structure match).
- P-44.4.1.8 / P-45.2: the -S(=O)NH2 locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as `_sulfinic_acid.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -S(=O)NH2 suffix. Confirmed via PubChem: 'cyclohexanesulfinamide'
  (O=S(N)C1CCCCC1).
- P-92/P-93.3.4.1 stereocenters: like `_sulfinic_acid.py`'s sulfur (and
  unlike `_sulfonamide.py`'s, whose two identical =O substituents keep it
  non-stereogenic), this module's sulfinamide sulfur (one =O, one N, one
  C, one lone pair -- four distinct "substituents") is itself a potential
  stereocenter in virtually every real -S(=O)NH2 molecule, confirmed via
  RDKit `FindPotentialStereo` on `CC(C)S(=O)N` (flags the sulfur even
  with no chain stereocenter at all). Same mechanism as
  `_sulfinic_acid.py`: `_common.heteroatom_stereo_prefix` cites a bare
  "(R)-"/"(S)-" prefix when the sulfur is the molecule's sole specified
  stereocenter (P-93.3.4.1); a specified sulfur stereocenter combined
  with a specified chain/ring carbon one remains out of scope.
- The sulfinamide nitrogen may carry zero, one, or two plain,
  unsubstituted, saturated, acyclic alkyl substituents (branched or
  unbranched), each cited as its own 'N-'-prefixed substituent directly
  ahead of the parent stem, in alphanumerical order (P-14.5.2, ignoring
  italicized prefixes like 'tert-' via `alpha_sort_key`), with a 'di'
  multiplying prefix (and a single shared 'N,N-' pair) when both are
  identical -- exactly `_sulfonamide.py`'s own N-/N,N-disubstitution
  pattern (same PR), built with `name_branch` (P-29 PIN style, fixed
  project-wide by PR #237). A *compound* N-substituent (has its own
  locant) is always parenthesized -- 'N-(propan-2-yl)methanesulfinamide',
  not PubChem's own raw 'N-propan-2-ylmethanesulfinamide' (CID 14896695),
  same correction as `_urea.py` (see that module's docstring for the
  Blue Book citations). Confirmed via PubChem for the non-compound case:
  'N-methylmethanesulfinamide' (CS(=O)NC).

Scope, deliberately narrow, mirroring `_sulfinic_acid.py`'s own first pass
exactly: a single -S(=O)NH2 on an acyclic chain or on a single saturated
carbon ring (monocyclic), with no other heteroatom anywhere in the molecule
except the sulfinamide group's own oxygen/nitrogen (and any N-alkyl
substituent's carbons) -- acid/amide-vs-other Table 3.3 seniority
coexistence is future work.

P-31.1.3: a monocyclic ring bearing a sulfinamide and exactly one C=C
ring double bond -- e.g. 'cyclohex-2-ene-1-sulfinamide' -- mirrors
`_sulfonamide.py`'s identical extension; the exact ring structures aren't
PubChem-registered, but the acyclic ene+sulfinamide combination is
('prop-2-ene-1-sulfinamide'), and the identical ring shape is confirmed
for the sibling `_sulfonic_acid.py`/`_sulfinic_acid.py`/`_sulfonamide.py`
modules. Deliberately narrow: a ring triple bond, any other ring
substituent alongside the ring double bond, and an N-substituent
alongside a ring double bond (an unverified combination) are all still
explicitly rejected pending further verification.

Explicitly out of scope (raise
`UnsupportedStructure`): an unsaturated or ring-bearing N-substituent (a
branched but otherwise plain saturated acyclic N-substituent is
supported, see above); an N-substituted sulfinamide, more than one
sulfinamide, or a specified stereocenter, on a von Baeyer polycyclic or
spiro skeleton (a single primary -S(=O)NH2 on such a skeleton is
supported, P-23/P-24 numbering integration via `_polycyclic_suffix.py`
with `elide_e=False`, mirroring `_thiol.py`); unsaturation reaching
outside the ring or a ring triple bond, a -S(=O)NH2 on a
substituent branch off an otherwise-unsubstituted *saturated* ring, two
or more -S(=O)NH2 groups, and a sulfinamide on a carbon that is also part
of a C=C/C#C bond. Two *aromatic*-ring cases:
`_name_phenyl_chain_sulfinamide` names a -S(=O)NH2 chain hanging off a
single plain, unsubstituted benzene ring (e.g.
'3-phenylpropane-1-sulfinamide'), mirroring `_sulfonamide.py`'s
identical benzene-ring-substituent path -- narrower than the acyclic
path: no N-alkyl substitution and no chain unsaturation.
`_name_benzenesulfinamide` names -S(=O)NH2 directly on a benzene ring
carbon (with or without other ring substituents), e.g.
'benzenesulfinamide' (PubChem PUG REST match for c1ccccc1S(=O)N),
mirroring `_sulfonamide.py`'s `_name_benzenesulfonamide` with the
retained name 'benzene' as stem -- the -S(=O)NH2's own locant is never
cited here; narrower than the acyclic path: no N-alkyl substitution
(a specified sole sulfur stereocenter is supported, see module docstring,
P-93.3.4.1).
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    all_chains,
    bfs,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    heteroatom_stereo_prefix,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    lowest_locant_set,
    most_multiple_bonds,
    name_from_substituents,
    non_single_bonds,
    ring_bond_locants,
    ring_chain_attachment,
    ring_branch_attachments,
    separate_aromatic_monocycles,
    ring_cycle,
    substituent_locant_set_and_citation,
    suffix_body,
)
from ._bicyclic import find_bicyclic_core
from ._numerals import alkane_name, alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._substituents import (
    format_substituent_prefixes,
    name_branch,
    substituents_for_ring,
    substituents_for_chain,
)

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _sulfinamide_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfinamide group: bonded to exactly one
    carbon, one double-bonded (terminal) oxygen, and one single-bonded
    primary amide nitrogen (terminal, two H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(oxygens) != 1 or len(nitrogens) != 1:
            continue
        (oxygen,) = oxygens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        if oxygen.GetDegree() != 1:
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


def has_sulfinamide_shape(mol) -> bool:
    return bool(_sulfinamide_sulfur_atoms(mol))


def _validate_and_collect_sulfinamides(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_sulfinamide`'s benzene-
    ring-substituent path (see `_name_phenyl_chain_sulfinamide`) can
    reuse this same validation for the rest of the molecule. Empty by
    default, so every other caller's behavior is unchanged."""
    sulfur_atoms = _sulfinamide_sulfur_atoms(mol)
    if not sulfur_atoms:
        raise UnsupportedStructure(
            "no sulfinamide (-S(=O)NH2) group found; this module only "
            "handles sulfinamides"
        )
    if len(sulfur_atoms) > 1:
        raise UnsupportedStructure(
            "more than one sulfinamide group is out of scope for this "
            "module"
        )
    sulfinamide_atom_idxs = set()
    for s in sulfur_atoms:
        sulfinamide_atom_idxs.add(s.GetIdx())
        sulfinamide_atom_idxs.update(
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
        elif atom.GetIdx() not in sulfinamide_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a sulfinamide group (P-65.3.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular an N-substituted sulfinamide, or a "
                "coexisting carboxylic/sulfonic/sulfinic acid or other "
                "characteristic group, needs handling not yet implemented "
                "here"
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


def _reject_enesulfinamide_carbon(graph, so_nh2_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if so_nh2_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a sulfinamide on a carbon that is also part of a C=C/C#C bond "
            "is out of scope for this module"
        )


def _n_alkyl_info(full_graph, full_carbon_graph, nitrogen_idx, n_alkyl_carbons, non_single_bond_atoms, halogens, mol=None):
    """Mirrors `_sulfonamide.py`'s own N-alkyl handling (PR #335 lineage):
    each N-substituent must be a plain, unsubstituted, saturated, acyclic
    alkyl chain (branched or unbranched). Returns (n_names,
    n_substituent_atoms) -- the latter must be excluded from the carbon
    graph before picking the principal chain, since these carbons hang
    off the (excluded) sulfinamide nitrogen rather than off the
    sulfinamide carbon itself."""
    n_names = []
    n_substituent_atoms = set()
    for n_alkyl_c in n_alkyl_carbons:
        n_atoms, _ = bfs(full_carbon_graph, n_alkyl_c)
        n_atoms = set(n_atoms)
        if n_atoms & non_single_bond_atoms:
            raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
        if any(nbr in n_atoms for h in halogens for nbr in full_graph[h]):
            raise UnsupportedStructure(
                "a substituted N-substituent (e.g. bearing a halogen) is "
                "not supported yet; only a plain, unsubstituted alkyl "
                "N-substituent is in scope"
            )
        n_names.append(name_branch(full_graph, n_alkyl_c, nitrogen_idx, {}, mol=mol))
        n_substituent_atoms |= n_atoms
    return n_names, n_substituent_atoms


def _name_from_substituents(chain_length, so_nh2_locant, ene_locants, yne_locants, grouped, n_names=()):
    # Only chain_length == 1 omits a substituent prefix's own locant too
    # (see `_alcohol.py`'s equivalent comment) -- an 'N-' locant is never
    # omitted either way (see `_add_n_names`/`format_substituent_prefixes`:
    # it marks a different atom than the chain itself).
    prefix = format_substituent_prefixes(_add_n_names(grouped, n_names), omit_locants=chain_length == 1)
    return prefix + name_from_substituents(chain_length, ene_locants, yne_locants, "sulfinamide", [so_nh2_locant], substituted=bool(grouped))


def _candidate_key(chain_length, so_nh2_locant, ene_locants, yne_locants, substituents, n_names=()):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, so_nh2_locant, ene_locants, yne_locants, grouped, n_names)
    return (
        (
            so_nh2_locant,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )

def _ring_name_from_substituents(ring_size, so_nh2_locant, ene_locants, yne_locants, grouped, n_names=()):
    has_unsaturation = bool(ene_locants or yne_locants)
    stem = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if not has_unsaturation:
        if total_subs == 0 and not n_names:
            # P-14.3.3: the sole substituent on an otherwise unsubstituted
            # ring has no locant to distinguish, e.g.
            # 'cyclohexanesulfinamide'.
            return stem + "sulfinamide"
        if total_subs == 0:
            prefix = format_substituent_prefixes(_add_n_names(grouped, n_names), omit_locants=True)
            return f"{prefix}{stem}sulfinamide"
        prefix = format_substituent_prefixes(_add_n_names(grouped, n_names))
        return f"{prefix}{stem}-{so_nh2_locant}-sulfinamide"

    # A competing ring double/triple bond (P-31.1.3) means the
    # sulfinamide's locant is never omittable even when it's the sole
    # substituent -- mirrors `_sulfonamide.py`'s identical treatment.
    unsaturated_stem = stem[:-3]
    needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    prefix = format_substituent_prefixes(_add_n_names(grouped, n_names))
    body = suffix_body(ene_locants, yne_locants, "sulfinamide", [so_nh2_locant])[0]
    return prefix + unsaturated_stem + ("a" if needs_stem_a else "") + "-" + body


def _ring_candidate_key(ring_size, so_nh2_locant, ene_locants, yne_locants, substituents, n_names=()):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, so_nh2_locant, ene_locants, yne_locants, grouped, n_names)
    return so_nh2_locant, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _name_cyclic_sulfinamide(mol, sulfur_idx, so_nh2_carbon, n_names, bonds=()):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    if bonds and any(substituents_for_ring(graph, ring_order, halogens, excluded, mol=mol).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "sulfinamide is not supported yet (see module docstring)"
        )
    if bonds and n_names:
        raise UnsupportedStructure(
            "an N-substituent alongside a ring double bond is not "
            "supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so_nh2_locant = position_of[so_nh2_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, so_nh2_locant, ene_locants, yne_locants, substituents, n_names)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _add_n_names(grouped, n_names):
    """Merge N-alkyl substituent names into a *display-only* copy of the
    ring `grouped` dict, each getting the non-numeric locant 'N' --
    mirrors `_amine.py`'s/`_sulfonamide.py`'s identical `_add_n_names`
    (P-66.4's N-locant interleaving, PR #443)."""
    if not n_names:
        return grouped
    display = {name: {"locants": list(info["locants"]), "compound": info["compound"]} for name, info in grouped.items()}
    for name, is_compound in n_names:
        info = display.setdefault(name, {"locants": [], "compound": is_compound})
        info["locants"].append("N")
    return display


def _benzenesulfinamide_name_from_substituents(grouped, n_names=()):
    # Mirrors `_sulfonamide.py`'s `_benzenesulfonamide_name_from_substituents`:
    # the mancude ring's own numbering is always free to start at the
    # -S(=O)NH2 carbon, so its locant is never cited even when other
    # substituents need theirs.
    display = _add_n_names(grouped, n_names)
    if not display:
        return "benzenesulfinamide"
    return f"{format_substituent_prefixes(display)}benzenesulfinamide"


def _benzenesulfinamide_candidate_key(so_nh2_locant, substituents, n_names=()):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzenesulfinamide_name_from_substituents(grouped, n_names)
    return so_nh2_locant, locant_set, citation_locants, name


def _name_benzenesulfinamide(mol, ring_atoms):
    """P-65.3.1: -S(=O)NH2 attached directly to a benzene ring carbon --
    e.g. 'benzenesulfinamide' (PubChem PUG REST match for
    c1ccccc1S(=O)N), and now 'N-methylbenzenesulfinamide' (PubChem
    PUG REST IUPACName match). Mirrors `_sulfonamide.py`'s
    `_name_benzenesulfonamide` with the retained name 'benzene' as stem;
    the -S(=O)NH2's own locant is never cited here. The N-alkyl
    substituent itself reuses `_n_alkyl_info` unchanged from the acyclic
    path (`name_sulfinamide`), merged into the ring citation via
    `_add_n_names` (so a coinciding name, e.g. 'N,4-dimethyl...',
    collapses into one multiplied citation like PubChem's own name). A
    specified sulfur stereocenter (with no ring-carbon one alongside it)
    gets a bare (R)-/(S)- prefix -- see module docstring, P-93.3.4.1."""
    sulfur_idx, so_nh2_carbon, nitrogen_idx, n_alkyl_carbons = _validate_and_collect_sulfinamides(
        mol, aromatic_ring_atoms=ring_atoms
    )
    stereo_prefix = heteroatom_stereo_prefix(mol, sulfur_idx) or ""

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
            so_nh2_locant = position_of[so_nh2_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _benzenesulfinamide_candidate_key(so_nh2_locant, substituents, n_names)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return stereo_prefix + best_name


def _name_phenyl_chain_sulfinamide(mol, ring_atoms):
    """Name a sulfinamide whose -S(=O)NH2 lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropane-1-sulfinamide.
    The ring is cited as a 'phenyl' substituent prefix on the chain,
    which is the parent hydride, mirroring `_sulfonamide.py`'s
    `_name_phenyl_chain_sulfonamide`. Narrower than the acyclic path
    above: no N-alkyl substitution and no chain unsaturation -- each is
    a separate follow-up (see
    tasks/phenyl-substituent-on-sulfinamide-chain.md's scope note)."""
    sulfur_idx, so_nh2_carbon, _, n_alkyl_carbons = _validate_and_collect_sulfinamides(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if n_alkyl_carbons:
        raise UnsupportedStructure(
            "an N-alkyl-substituted sulfinamide alongside a benzene-ring "
            "substituent is not supported yet"
        )
    stereo_prefix = heteroatom_stereo_prefix(mol, sulfur_idx) or ""
    excluded = {sulfur_idx}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if sulfur_idx not in (b[0], b[1]) and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "sulfinamide chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_branch_attachments(mol, graph, rings)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain sulfinamide is not "
            "supported yet"
        )
    chain, branches = longest_branched_chain_through(graph, so_nh2_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        so_nh2_locant = position_of[so_nh2_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, so_nh2_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return stereo_prefix + best_name


def _name_ring_substituent_chain_sulfinamide(mol, sulfur_idx, so_nh2_carbon):
    """Name a sulfinamide whose -S(=O)NH2 lies entirely on a single
    branched chain hanging off one atom of an otherwise-plain saturated
    monocyclic ring (the ring itself bears no sulfinamide) -- e.g.
    cyclohexylmethanesulfinamide. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_sulfinamide` above and
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

    chain, branches = longest_branched_chain_through(graph, so_nh2_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
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
        so_nh2_locant = position_of[so_nh2_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, so_nh2_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_von_baeyer_or_spiro_sulfinamide(mol, sulfur_idx, so_nh2_carbon, n_alkyl_carbons, bonds, stereo_prefix):
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single sulfinamide (-S(=O)NH2) suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`'s
    `elide_e=False` ('sulfinamide' begins with a consonant, P-16.3.3, same
    as `_thiol.py`'s 'thiol' -- e.g. 'bicyclo[2.2.1]heptane-2-
    sulfinamide', PubChem CID 165459501). Mirrors `_thiol.py`'s own
    bicyclic/polycyclic-before-spiro dispatch order and restrictions:
    exactly one primary sulfinamide on the ring system itself, no
    specified stereocenter (on the ring or the sulfinamide sulfur
    itself); ring unsaturation composes on the bicyclic/polycyclic
    branch (matches `_ketone.py`'s bicyclic/polycyclic-vs-monospiro split) but not the
    monospiro branch, which has no base mechanism yet (P-31.1.5)."""
    if n_alkyl_carbons:
        raise UnsupportedStructure(
            "an N-substituted sulfinamide on a von Baeyer bicyclic/"
            "polycyclic or monospiro ring system is not supported yet"
        )
    if stereo_prefix:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a von Baeyer bicyclic/"
            "polycyclic or monospiro sulfinamide is not supported yet "
            "(see P-92)"
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
            so_nh2_carbon,
            {sulfur_idx},
            "sulfinamide",
            "sulfinamide",
            bicyclic_core,
            polycyclic_core,
            von_baeyer_ring_count,
            elide_e=False,
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
            mol, so_nh2_carbon, {sulfur_idx}, "sulfinamide", "sulfinamide", spiro_atom, elide_e=False
        )

    raise UnsupportedStructure(
        "polycyclic and fused-ring sulfinamides are not supported yet "
        "(P-23/P-25 numbering integration with a suffix group is future "
        "work)"
    )


def name_sulfinamide(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        return _name_phenyl_chain_sulfinamide(mol, set().union(*aromatic_rings))
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            _, so_nh2_carbon, _, _ = _validate_and_collect_sulfinamides(mol, aromatic_ring_atoms=ring_atoms)
            if so_nh2_carbon in ring_atoms:
                return _name_benzenesulfinamide(mol, ring_atoms)
            return _name_phenyl_chain_sulfinamide(mol, ring_atoms)
    sulfur_idx, so_nh2_carbon, nitrogen_idx, n_alkyl_carbons = _validate_and_collect_sulfinamides(mol)
    stereo_prefix = heteroatom_stereo_prefix(mol, sulfur_idx) or ""
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
    _reject_enesulfinamide_carbon(graph, so_nh2_carbon, bonds)
    non_single_bond_atoms = {a for a, b, _ in all_non_single} | {b for a, b, _ in all_non_single}
    halogens = halogen_substituents(mol)
    n_names, n_substituent_atoms = _n_alkyl_info(
        graph, carbon_adjacency(mol), nitrogen_idx, n_alkyl_carbons, non_single_bond_atoms, halogens, mol=mol
    )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        return _name_von_baeyer_or_spiro_sulfinamide(
            mol, sulfur_idx, so_nh2_carbon, n_alkyl_carbons, bonds, stereo_prefix
        )
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
        if any(order == _YNE_ORDER for _, _, order in ring_bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a sulfinamide "
                "is not supported yet -- only a ring double bond is in "
                "scope for this first pass (see P-31.1.3)"
            )
        if so_nh2_carbon not in ring_atoms:
            if not bonds and not n_alkyl_carbons:
                return stereo_prefix + _name_ring_substituent_chain_sulfinamide(mol, sulfur_idx, so_nh2_carbon)
            raise UnsupportedStructure(
                "a sulfinamide on a substituent branch chain rather than "
                "the ring itself is not supported yet"
            )
        return stereo_prefix + _name_cyclic_sulfinamide(mol, sulfur_idx, so_nh2_carbon, n_names, bonds)

    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    full_carbon_graph = carbon_adjacency(mol)
    carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in n_substituent_atoms}
    chains = all_chains(carbon_graph)

    eligible = []
    for chain in chains:
        if so_nh2_carbon not in chain:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the sulfinamide-bearing carbon (and/or a multiple bond) does "
            "not lie on a single longest carbon chain; a shorter principal "
            "chain is not supported yet"
        )

    best_key = None
    best_name = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so_nh2_locant = position_of[so_nh2_carbon]
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, so_nh2_locant, ene_locants, yne_locants, substituents, n_names)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return stereo_prefix + best_name
