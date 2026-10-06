"""Naming of nitriles (the '-nitrile' suffix, terminal -C#N) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-66.5, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  'nitrile' is the preselected suffix for -C#N, cited in a combined
  chain/suffix name the same way 'al'/'one'/'ol'/'amine' are for their own
  groups (`_aldehyde.py`/`_ketone.py`/`_alcohol.py`/`_amine.py`), e.g.
  'propanenitrile', 'pent-4-enenitrile'. Table 3.3 ranks 'nitrile' junior to
  'amide' and senior to 'al'; this module never attempts suffix-vs-suffix
  seniority competition against a coexisting senior or junior group, since
  only C, N (the nitrile nitrogen itself), and halogen atoms are accepted at
  all (any oxygen routes elsewhere in `core.py`, or is out of scope here).
- Unlike 'al'/'one'/'amine', 'nitrile' begins with a consonant, so the
  euphonic elision `_aldehyde.py` applies before a vowel-initial suffix
  never triggers here: the stem's final 'e' (from the parent alkane/alkene/
  alkyne name) is always kept, e.g. 'propanenitrile', not 'propannitrile';
  'pent-4-enenitrile', not 'pent-4-ennitrile'.
- Like -CHO (`_aldehyde.py`), a -C#N carbon is always a chain terminus (it
  has exactly one carbon neighbor, being otherwise saturated by the triple
  bond to N), so it is never a genuine locant choice: whichever end of the
  principal chain bears it is numbered C1 (P-44.4.1.8, suffix locants
  minimized first), and that locant is always '1' and never cited (P-14.3.3):
  'propanenitrile', not 'propane-1-nitrile'.
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-nitrile name reuses the same mechanics as
  `_unsaturated.py`/`_aldehyde.py` (ending replaces 'ane' entirely,
  multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before
  'yne', euphonic stem 'a' before a multiplied ending); the 'nitrile' suffix
  word is appended directly onto the last such segment, since it never
  carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'nitrile' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than two nitrile groups - a third can't sit on both chain termini
  the way exactly two can (P-66.6.3's dinitrile case, e.g.
  'butanedinitrile', is supported for the acyclic case), and this project
  has no 'cyano' substituent-prefix support yet for a branch-mounted one.
  A dinitrile alongside any ring is also still out of scope.
- Any oxygen - routed to a different module by `core.py`, or out of scope
  entirely if this module is called directly on one.
- A nitrile nitrogen that isn't a plain, isolated -C#N (any degree other
  than 1, or a bond order other than 3.0 to its one carbon neighbor).
- A nitrile carbon with other than exactly one carbon neighbor: zero (bare
  HC#N) or two-or-more (not a valid nitrile shape) - neither is a terminal
  substitutive -C#N.
- An aromatic nitrile carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory, except for one narrow case: a
  nitrile's chain hanging off a single plain, unsubstituted benzene ring
  with no other substituent on the ring (`_name_phenyl_chain_nitrile`,
  e.g. '3-phenylpropanenitrile'), mirroring `_aldehyde.py`/`_amide.py`'s
  identical benzene-ring-substituent path. Narrower than the acyclic path:
  no chain unsaturation, no specified stereocenter, and no substituted
  benzene/naphthalene - each a separate follow-up.
- A ring shape other than a single saturated monocyclic all-carbon ring
  or a single plain benzene ring, with -C#N directly on it -- P-66.5.1.1.3's
  'carbonitrile' suffix (a substituent-style name rather than this
  module's own parent-chain suffix) is handled by `_name_ring_nitrile`/
  `_name_benzonitrile` for those two single-ring cases only.
- Any other heteroatom (O, S, ...), or a nitrile carbon entangled with
  another heteroatom.

P-91.3/P-92: a molecule with one
or more *specified* tetrahedral stereocenters -- every one on the
principal chain itself, no unspecified one alongside them, and no
C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
ascending locant order, same pattern as `_amine.py`/`_thiol.py` (chain
only, this module has no ring support). The nitrile carbon itself (sp,
triple-bonded to nitrogen) is never a potential stereocenter, confirmed
via RDKit `FindPotentialStereo` on `CC[C@@H](C)C#N`.
"""

from rdkit import Chem

from ._common import (
    stereo_locant_rank,
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    all_chains,
    carbon_adjacency,
    carbon_on_ring,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    heteroaromatic_monocycle_name,
    is_plain_benzene_ring,
    longest_branched_chain,
    lowest_locant_set,
    most_multiple_bonds,
    multiplied_word,
    name_from_substituents,
    non_single_bonds,
    ring_branch_attachments,
    ring_hosting_anchors,
    separate_aromatic_monocycles,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereo_elements,
    specified_stereocenters,
    stereo_element_atoms,
    substituent_locant_set_and_citation,
    two_separate_rings_with_plain_aromatic_substituent,
)
from ._substituents import (
    substituents_for_chain,
    format_substituent_prefixes,
    name_branch,
)

_NITRILE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def has_nitrile_shape(mol) -> bool:
    """True iff `mol` has at least one nitrogen atom of degree 1, triple-
    bonded to a carbon (a plain -C#N pattern), regardless of whether the
    rest of the molecule is in scope. Used by `core.py` to route ahead of
    the amine dispatch within the no-oxygen nitrogen branch, since a
    nitrile nitrogen would otherwise look like an unhandled shape to
    `name_amine`'s own bond-order check."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        (neighbor,) = atom.GetNeighbors()
        if bond.GetBondTypeAsDouble() == _NITRILE_ORDER and neighbor.GetAtomicNum() == 6:
            return True
    return False


def _validate_and_collect_nitriles(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return the set of nitrile-nitrogen atom indices.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring or heteroaromatic monocycle (pyridine/furan/thiophene/pyrrole)
    with exactly one exocyclic attachment -- exempted wholesale from the
    per-atomic-number checks below (already independently verified by that
    shape check itself) so `name_nitrile`'s benzene/heteroaromatic-ring-
    substituent path (see `_name_phenyl_chain_nitrile`) and its
    `two_separate_rings_with_plain_aromatic_substituent` shape (mirroring
    `_thiol.py`'s #628) can both reuse this same validation for the rest
    of the molecule. Empty by default, so every other caller's behavior is
    unchanged."""
    nitriles = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atom.GetIdx() in aromatic_ring_atoms:
            if atomic_num == 6:
                has_carbon = True
            continue
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a nitrile nitrogen (P-66.5) and "
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
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a nitrogen bonded to more than one heavy atom is not a "
                    "plain nitrile (-C#N) and is out of scope for this "
                    "module (P-66.5)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a nitrile nitrogen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() != _NITRILE_ORDER:
                raise UnsupportedStructure(
                    "a nitrogen not triple-bonded to its one carbon neighbor "
                    "is not a nitrile and is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a nitrile carbon on an aromatic ring is out of scope "
                    "for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                raise UnsupportedStructure(
                    "a nitrile carbon with other than exactly one carbon "
                    "neighbor is out of scope for this module (Table 3.3)"
                )
            nitriles.add(atom.GetIdx())
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
    if not nitriles:
        raise UnsupportedStructure("no nitrile (-C#N) group found; this module only handles nitriles")
    if len(nitriles) > 2:
        raise UnsupportedStructure(
            "more than two nitrile groups is not supported yet (a third "
            "nitrile can't sit on both chain termini, and this project "
            "has no 'cyano' substituent-prefix support yet)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return nitriles


def _name_from_substituents(chain_length, nitrile_count, ene_locants, yne_locants, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(
        chain_length, ene_locants, yne_locants, multiplied_word(nitrile_count, "nitrile")
    )


def _candidate_key(chain_length, nitrile_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first. The nitrile locant set still drives
    orientation choice even though it is never printed (see module
    docstring)."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    nitrile_locant_set = lowest_locant_set(nitrile_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, len(nitrile_locants), ene_locants, yne_locants, grouped)
    return (
        (
            nitrile_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _nitrile_locants(position_of, nitriles, graph):
    locants = []
    for n in nitriles:
        (carbon,) = graph[n]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _name_acyclic_nitrile(mol, nitriles, bonds, stereo=None):
    """`stereo`: None, or a list of ("atom"/"bond", idx, "R"/"S"/"E"/"Z")
    from `specified_stereo_elements`/`specified_stereocenters` -- if given,
    only chain candidates that include every tetrahedral stereocenter are
    eligible (P-92: a stereocenter on a substituent branch is out of
    scope; a double-bond E/Z element's atoms are already required to lie
    on the chain via `bonds`, so no separate check is needed for those),
    and the winning candidate's own locants are used to format a
    "(<locant><R/S/E/Z>,...)-" prefix onto the final name, ascending
    locant order (P-91.3, including when both kinds coexist)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = all_chains(carbon_adjacency(mol))
    stereo_atoms = stereo_element_atoms(mol, stereo) if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _nitrile_locants(position_of, nitriles, graph) is None:
            continue
        if stereo is not None and any(atom not in position_of for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the nitrile-bearing carbon (and/or multiple bond, and/or a "
            "stereocenter on a substituent branch, see P-92) does not lie "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing the -C#N group (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            nitrile_locants = _nitrile_locants(position_of, nitriles, graph)
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, nitriles, mol=mol)
            key, name = _candidate_key(chain_length, nitrile_locants, ene_locants, yne_locants, substituents)
            key = (key, stereo_locant_rank(mol, stereo, position_of))
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = []
        for kind, idx, code in stereo:
            if kind == "atom":
                locant = best_position_of[idx]
            else:
                bond = mol.GetBondWithIdx(idx)
                locant = min(best_position_of[bond.GetBeginAtomIdx()], best_position_of[bond.GetEndAtomIdx()])
            labels.append((locant, code))
        labels.sort()
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _ring_substituents(graph, ring_order, halogens, excluded, mol=None, aromatic_atoms=frozenset()):
    """{ring position -> [substituent name, ...]}, mirroring
    `_aldehyde.py`'s identically-named helper -- every branch hanging off
    a ring atom other than the nitrile carbon itself (in `excluded`) is a
    plain substituent prefix (alkyl/halogen), or the aromatic ring itself
    (phenyl or a heteroaromatic monocycle, via `name_branch`) when reached
    as one half of a `two_separate_rings_with_plain_aromatic_substituent`
    shape."""
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [
                name_branch(graph, root, atom, halogens, aromatic_atoms, mol=mol) for root in branch_roots
            ]
    return substituents


def _ring_name_from_substituents(ring_size, cn_locant, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(ring_size, [], [], prefix, total_subs, "carbonitrile", [cn_locant])


def _ring_candidate_key(ring_size, cn_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _ring_name_from_substituents(ring_size, cn_locant, grouped)
    return cn_locant, locant_set, citation_locants, name


def _name_ring_nitrile(mol, ring_atoms, stereo=None, aromatic_atoms=frozenset()):
    """P-66.5.1.1.3: "The suffix 'carbonitrile' is always used to name
    nitriles having the -C#N group attached to a ring" -- e.g.
    'cyclohexanecarbonitrile (PIN)'. Unlike the chain-parent 'nitrile'
    suffix above (where the -C#N carbon is always the parent's own C1),
    here the ring itself is the parent hydride and the -C#N carbon is a
    substituent atom hanging directly off one ring carbon, mirroring
    `_aldehyde.py`'s `_name_ring_aldehyde` construction.

    A saturated monocyclic all-carbon ring with exactly one -C#N hanging
    directly off one ring atom (no ring unsaturation), plus any number of
    substituents (alkyl/halogen) on any ring atom, including the -C#N-
    bearing atom itself.

    `stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- every stereocenter must lie on the ring
    itself (a stereocenter on the -C#N substituent branch is out of
    scope for now).

    `aromatic_atoms`: when the nitrile-bearing ring reaches this function
    as one half of a `two_separate_rings_with_plain_aromatic_substituent`
    shape (mirroring `_thiol.py`'s #628), the other, aromatic ring's atoms
    -- cited as a plain substituent (phenyl or a heteroaromatic monocycle)
    via `name_branch` the same way an ordinary alkyl ring substituent
    already is. Empty by default, so the original single-ring dispatch is
    unchanged."""
    (nitrile_nitrogen,) = _validate_and_collect_nitriles(mol, aromatic_ring_atoms=aromatic_atoms)

    graph = adjacency(mol)
    (nitrile_carbon,) = graph[nitrile_nitrogen]
    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] != nitrile_nitrogen
        and b[1] != nitrile_nitrogen
        and b[0] in ring_atoms
        and b[1] in ring_atoms
        and not (b[0] in aromatic_atoms and b[1] in aromatic_atoms)
    ]
    if all_non_single:
        raise UnsupportedStructure(
            "an unsaturated ring alongside a nitrile substituent is not "
            "supported yet (see P-31.1.3)"
        )

    ring_neighbors = [n for n in graph[nitrile_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "a nitrile not directly attached to a single ring atom is not "
            "supported yet"
        )
    (ring_atom,) = ring_neighbors

    if stereo is not None and any(atom not in ring_atoms for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on the -C#N substituent branch rather than "
            "the ring itself is not supported yet (see P-92)"
        )

    halogens = halogen_substituents(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            cn_locant = position_of[ring_atom]
            substituents = _ring_substituents(
                graph, candidate, halogens, {nitrile_carbon}, mol=mol, aromatic_atoms=aromatic_atoms
            )
            key = _ring_candidate_key(ring_size, cn_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _benzonitrile_name_from_substituents(grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if total_subs == 0:
        # 'benzonitrile' is a fully retained name (P-66.5.1.1.3) -- unlike
        # 'cyclohexanecarbonitrile', there is no locant position to even
        # omit.
        return "benzonitrile"
    return f"{format_substituent_prefixes(grouped)}benzonitrile"


def _benzonitrile_candidate_key(cn_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzonitrile_name_from_substituents(grouped)
    return cn_locant, locant_set, citation_locants, name


def _name_benzonitrile(mol, ring_atoms, exempt_atoms=None):
    """P-66.5.1.1.3: 'benzonitrile' is a retained name that is itself the
    PIN for a -C#N hanging directly off one carbon of an otherwise-plain
    (or substituted) benzene ring -- e.g. 'benzonitrile' (PubChem CID
    7505), '2-methylbenzonitrile' (CID 10287), '4-methylbenzonitrile'
    (CID 12174). Structurally identical to `_name_ring_nitrile`'s
    saturated-ring case, mirroring `_aldehyde.py`'s `_name_benzaldehyde`,
    with the retained name 'benzonitrile' replacing 'cyclo' + alkane_name
    + 'carbonitrile' as the whole suffix unit (no locant is ever cited
    for the -C#N position itself)."""
    (nitrile_nitrogen,) = _validate_and_collect_nitriles(mol, aromatic_ring_atoms=exempt_atoms or ring_atoms)

    graph = adjacency(mol)
    (nitrile_carbon,) = graph[nitrile_nitrogen]
    ring_neighbors = [n for n in graph[nitrile_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "a nitrile not directly attached to a single ring atom is not "
            "supported yet"
        )
    (ring_atom,) = ring_neighbors
    other_ring_atom_branches = [n for n in graph[ring_atom] if n not in ring_atoms and n != nitrile_carbon]
    if other_ring_atom_branches:
        raise UnsupportedStructure(
            "a substituent on the same ring atom as the nitrile is not "
            "supported yet"
        )

    halogens = halogen_substituents(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            cn_locant = position_of[ring_atom]
            substituents = _ring_substituents(graph, candidate, halogens, {nitrile_carbon}, mol=mol)
            key = _benzonitrile_candidate_key(cn_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    return best_name


def _name_phenyl_chain_nitrile(mol, ring_atoms):
    """Name a nitrile whose -C#N lies entirely on a single unbranched chain
    hanging off one atom of an otherwise-plain, unsubstituted benzene ring
    -- e.g. 3-phenylpropanenitrile. The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring
    `_aldehyde.py`'s `_name_phenyl_chain_aldehyde`. Narrower than the
    acyclic path above: no chain unsaturation and no specified stereocenter
    -- each is a separate follow-up."""
    nitriles = _validate_and_collect_nitriles(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "nitrile chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in nitriles and b[1] not in nitriles and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "nitrile chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_branch_attachments(mol, graph, rings)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain nitrile is not "
            "supported yet"
        )
    (nitrile_nitrogen,) = nitriles
    (nitrile_carbon,) = graph[nitrile_nitrogen]
    chain, branches = longest_branched_chain(graph, nitrile_carbon, ring_atoms, nitriles, halogens=halogen_substituents(mol))

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, 1, [], [], grouped)


def name_nitrile(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        union = set().union(*aromatic_rings)
        anchors = [
            next(iter(adjacency(mol)[n])) for n in _validate_and_collect_nitriles(mol, aromatic_ring_atoms=union)
        ]
        host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
        if host is not None:
            return _name_benzonitrile(mol, host, union)
        return _name_phenyl_chain_nitrile(mol, union)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        is_benzene = is_plain_benzene_ring(mol, ring_atoms)
        is_heteroaromatic = not is_benzene and (
            heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None
        )
        if is_benzene or is_heteroaromatic:
            ring_nitriles = _validate_and_collect_nitriles(mol, aromatic_ring_atoms=ring_atoms)
            if len(ring_nitriles) != 1:
                raise UnsupportedStructure(
                    "more than one nitrile group alongside a ring is not "
                    "supported yet (only the acyclic dinitrile case is)"
                )
            (only_nitrogen,) = ring_nitriles
            graph = adjacency(mol)
            (only_carbon,) = graph[only_nitrogen]
            ring_neighbors = [n for n in graph[only_carbon] if n in ring_atoms]
            if len(ring_neighbors) == 1:
                if is_benzene:
                    return _name_benzonitrile(mol, ring_atoms)
                raise UnsupportedStructure(
                    "a nitrile directly attached to a heteroaromatic ring "
                    "is not supported yet (only a plain benzene ring's "
                    "'benzonitrile' direct-attachment path, or a "
                    "heteroaromatic ring cited as a chain substituent, is)"
                )
            return _name_phenyl_chain_nitrile(mol, ring_atoms)
        ring_nitriles = _validate_and_collect_nitriles(mol)
        if len(ring_nitriles) != 1:
            raise UnsupportedStructure(
                "more than one nitrile group alongside a ring is not "
                "supported yet (only the acyclic dinitrile case is)"
            )
        (only_nitrogen,) = ring_nitriles
        graph = adjacency(mol)
        (only_carbon,) = graph[only_nitrogen]
        ring_neighbors = [n for n in graph[only_carbon] if n in ring_atoms]
        if len(ring_neighbors) == 1:
            return _name_ring_nitrile(mol, ring_atoms, specified_stereocenters(mol))

    # Two separate simple monocycles joined by one direct bond, one a plain
    # benzo/heteroaromatic ring with no substituent of its own (mirroring
    # `_thiol.py`'s #628) -- e.g. 2-phenylcyclohexane-1-carbonitrile --
    # reuses `_name_ring_nitrile` with the aromatic ring's atoms passed
    # through as an exemption/substituent. Narrower than that: only the
    # "every nitrile is ring-borne" split is handled below (matching this
    # shape's only verified real structures), mirroring #622/#624/#628's
    # identical scoping decision.
    if ring_info.NumRings() == 2:
        aromatic_shape = two_separate_rings_with_plain_aromatic_substituent(mol, adjacency(mol))
        if aromatic_shape is not None:
            ring_atoms, aromatic_atoms, _, _ = aromatic_shape
            ring_nitriles = _validate_and_collect_nitriles(mol, aromatic_ring_atoms=aromatic_atoms)
            if len(ring_nitriles) != 1:
                raise UnsupportedStructure(
                    "more than one nitrile group alongside this two-ring "
                    "aromatic-substituent shape is not supported yet"
                )
            (only_nitrogen,) = ring_nitriles
            graph = adjacency(mol)
            (only_carbon,) = graph[only_nitrogen]
            ring_neighbors = [n for n in graph[only_carbon] if n in ring_atoms]
            if len(ring_neighbors) != 1:
                raise UnsupportedStructure(
                    "a nitrile on a chain hanging off the non-aromatic "
                    "ring, with the ring itself bearing no nitrile of its "
                    "own, alongside this two-ring aromatic-substituent "
                    "shape is not supported yet"
                )
            return _name_ring_nitrile(
                mol, ring_atoms, specified_stereocenters(mol), aromatic_atoms=aromatic_atoms
            )

    nitriles = _validate_and_collect_nitriles(mol)
    stereo = specified_stereo_elements(mol)
    graph = adjacency(mol)
    # Exclude each C#N nitrile bond itself: `non_single_bonds` reports it as
    # order 3.0 same as a C#C, but it isn't a chain 'yne' bond (one endpoint
    # is the nitrile nitrogen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in nitriles and b[1] not in nitriles]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    if carbon_on_ring(mol, adjacency(mol), [c for n in nitriles for c in adjacency(mol)[n]]):
        raise UnsupportedStructure(
            "this ring shape alongside a nitrile (ring unsaturation, or a "
            "ring other than a single saturated monocyclic/benzene one) is "
            "out of scope for this module's 'carbonitrile' suffix path "
            "(P-66.5.1.1.3)"
        )
    return _name_acyclic_nitrile(mol, nitriles, bonds, stereo)
