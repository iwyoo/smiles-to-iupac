"""Naming of aldehydes (the '-al' suffix, terminal -CHO) on acyclic saturated
or unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-33.3, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'al' is the preselected suffix for -CHO, cited in a combined chain/suffix
  name the same way 'ol'/'amine'/'one' are for -OH/-NH2/=O (`_alcohol.py`/
  `_amine.py`/`_ketone.py`), e.g. 'propanal', 'pentanedial'. Table 3.3 ranks
  'al' senior to 'one'/'ol'/'amine'; a carbonyl carbon shaped like a ketone
  (two carbon neighbors) or a carboxylic acid/ester/amide (a second oxygen on
  the same carbon) is rejected outright rather than silently named as if it
  were an aldehyde, and this module never attempts suffix-vs-suffix
  seniority competition against a coexisting senior group either (carboxylic
  acid/ester/amide/nitrile) since only C, halogen, and aldehyde/hydroxyl
  oxygen atoms are accepted at all. A coexisting hydroxyl (-OH), being junior
  to 'al', is *not* rejected: it is cited as the 'hydroxy' substituent prefix
  instead (P-41), e.g. 'OCC=O' -> '2-hydroxyethanal'.
- Unlike -OH/=O, a -CHO carbon is always a chain terminus (it has exactly one
  carbon neighbor, being otherwise saturated by =O and one H), so it is
  never a genuine locant choice: whichever end of the principal chain bears
  it is numbered C1 (P-44.4.1.8, suffix locants minimized first, same
  ordering as `_alcohol.py`/`_amine.py`/`_ketone.py`), and that locant is
  always '1' (or '1' and the last position, for a chain with -CHO at both
  ends) - so unlike those modules, this one never actually cites an 'al'
  locant in the name (P-14.3.3's general omission-when-unambiguous applies
  uniformly here, not just to a mononuclear/two-carbon parent): 'propanal',
  not 'propan-1-al'; 'pentanedial', not 'pentane-1,5-dial'. An 'ene'/'yne'
  locant elsewhere on the chain is still cited as usual, e.g. 'hex-4-enal'.
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-aldehyde name reuses the same mechanics
  as `_unsaturated.py`/`_ketone.py` (ending replaces 'ane' entirely,
  multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before
  'yne', euphonic stem 'a' before a multiplied ending); the 'al'/'dial'
  suffix word is appended directly onto the last such segment (eliding a
  trailing vowel the same way 'ene'+'ol' does in `_alcohol.py`), since it
  never carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'al' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.
- P-91.3/P-92: a molecule with one
  or more *specified* tetrahedral stereocenters -- every one on the
  principal chain itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
  ascending locant order, e.g. '(2R)-2-chloropropanal',
  '(2R,3S)-2,3-dichlorobutanal', same pattern as `_carboxylic_acid.py`
  (CIP computation delegated entirely to RDKit).

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that isn't a doubly-bonded, isolated aldehyde carbonyl oxygen or
  a singly-bonded hydroxyl (ethers, and any oxygen bonded to more than one
  heavy atom).
- A carbonyl carbon with other than exactly one carbon neighbor: zero (a
  carbon-less carbonyl, e.g. formaldehyde) or two-or-more (a ketone) -
  neither is this module's territory.
- An aromatic carbonyl carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory.
- Any other heteroatom (N, S, ...).
- A hydroxyl on a carbon that is also part of a C=C/C#C bond (an enol,
  tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
- More than one -CHO on a ring, a ring with other unsaturation, a
  standalone hydroxyl alongside a ring aldehyde, or any ring shape other
  than a single saturated monocyclic all-carbon ring or a single plain
  benzene ring -- P-66.6.1.1.3's 'carbaldehyde' suffix (a substituent-
  style name rather than this module's own parent-chain suffix) is
  handled by `_name_ring_aldehyde`/`_name_benzaldehyde` for the single-
  group case only.
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
    bond_locant,
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
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    ring_chain_attachments_with_halogens,
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
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
    substituents_for_chain,
)

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _validate_and_collect_aldehydes(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (aldehydes, hydroxyls): the set of aldehyde carbonyl-oxygen
    atom indices, and the set of any coexisting hydroxyl-oxygen atom indices.
    A hydroxyl is junior to 'al' in Table 3.3's suffix seniority order, so it
    is cited as the 'hydroxy' substituent prefix instead of competing for the
    suffix (P-41).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring or heteroaromatic monocycle (pyridine/furan/thiophene/pyrrole)
    with exactly one exocyclic attachment -- exempted wholesale from the
    per-atomic-number checks below (already independently verified by that
    shape check itself) so `name_aldehyde`'s benzene/heteroaromatic-ring-
    substituent path (see `_name_phenyl_chain_aldehyde`) and its
    `two_separate_rings_with_plain_aromatic_substituent` shape (mirroring
    `_nitrile.py`'s #631) can both reuse this same validation for the rest
    of the molecule. Empty by default, so every other caller's behavior is
    unchanged."""
    aldehydes = set()
    hydroxyls = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atom.GetIdx() in aromatic_ring_atoms:
            if atomic_num == 6:
                has_carbon = True
            continue
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an aldehyde carbonyl oxygen (P-33.3) "
                "and halogen substituents (P-35.2.1) are not supported yet"
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
                    "ether) is out of scope; only an isolated aldehyde "
                    "carbonyl or hydroxyl is supported (Table 3.3, P-33.3)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("an aldehyde/hydroxyl oxygen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() == 1.0:
                if atom.GetTotalNumHs() != 1:
                    raise UnsupportedStructure(
                        "an oxygen that isn't a carbonyl (=O) or hydroxyl "
                        "(-OH) is out of scope for this module"
                    )
                hydroxyls.add(atom.GetIdx())
                continue
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring is out of scope for this "
                    "module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                raise UnsupportedStructure(
                    "a carbonyl carbon with other than exactly one carbon "
                    "neighbor (a ketone, or a carbon-less carbonyl such as "
                    "formaldehyde) is not an aldehyde and is out of scope "
                    "for this module (Table 3.3)"
                )
            if mol.GetBondBetweenAtoms(carbon.GetIdx(), carbon_neighbors[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "a carbonyl carbon that is itself doubly bonded to its "
                    "carbon neighbor (a cumulated double bond, e.g. a "
                    "ketene's C=C=O) is not an aldehyde and is out of scope "
                    "for this module"
                )
            aldehydes.add(atom.GetIdx())
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
    if not aldehydes:
        raise UnsupportedStructure("no aldehyde (-CHO) group found; this module only handles aldehydes")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return aldehydes, hydroxyls


def _name_from_substituents(chain_length, al_count, ene_locants, yne_locants, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(
        chain_length, ene_locants, yne_locants, multiplied_word(al_count, "al")
    )


def _candidate_key(chain_length, al_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first. The 'al' locant set still drives
    orientation choice even though it is never printed (see module
    docstring)."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    al_locant_set = lowest_locant_set(al_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, len(al_locants), ene_locants, yne_locants, grouped)
    return (
        (
            al_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _al_locants(position_of, aldehydes, graph):
    locants = []
    for o in aldehydes:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants




def _name_acyclic_aldehyde(
    mol, aldehydes, hydroxyls, bonds, stereo=None, extra_names=None, required_atoms=frozenset(), carbon_graph=None
):
    """`stereo`: None, or a list of ("atom"/"bond", idx, "R"/"S"/"E"/"Z")
    from `specified_stereo_elements`/`specified_stereocenters` -- if given,
    only chain candidates that include every tetrahedral stereocenter are
    eligible (P-92: a stereocenter on a substituent branch rather than the
    principal chain is out of scope, mirroring `_carboxylic_acid.py`'s
    identical treatment; a double-bond E/Z element's atoms are already
    required to lie on the chain via `bonds`, so no separate check is
    needed for those), and the winning candidate's own locants are used to
    format a "(<locant><R/S/E/Z>,...)-" prefix onto the name, ascending
    locant order (P-91.3, including when both kinds coexist).

    `extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group demoted to a substituent prefix by
    `_seniority.senior_class` (e.g. a demoted amine's 'amino'), reused by
    `_coexisting_groups.py` so a pairwise module doesn't have to
    reimplement this function's chain search/candidate selection.
    `required_atoms`: additional carbon atoms (e.g. every demoted amine's
    own carbon neighbor) that a candidate chain must also carry -- both
    empty/None by default so existing callers are unaffected.
    `carbon_graph`: the carbon-only graph to search for the principal
    chain -- defaults to `carbon_adjacency(mol)` (unchanged behavior);
    `_ether_aldehyde.py` passes one with a coexisting ether's alkoxy-
    branch component already removed, the same reason and pattern as
    `_thiol.py`/`_ketone.py`'s identical `carbon_graph` parameter (PR
    #428/#429)."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}, **(extra_names or {})}
    chains = all_chains(carbon_graph if carbon_graph is not None else carbon_adjacency(mol))
    stereo_atoms = stereo_element_atoms(mol, stereo) if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _al_locants(position_of, aldehydes, graph) is None:
            continue
        if not required_atoms <= set(chain):
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _al_locants({a: i + 1 for i, a in enumerate(c)}, aldehydes, graph) is not None
            and required_atoms <= set(c)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereo element on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every aldehyde-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more -CHO groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            al_locants = _al_locants(position_of, aldehydes, graph)
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, aldehydes, mol=mol)
            key, name = _candidate_key(chain_length, al_locants, ene_locants, yne_locants, substituents)
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
    `_carboxylic_acid.py`'s identically-named helper -- every branch
    hanging off a ring atom other than the aldehyde carbon itself (in
    `excluded`) is a plain substituent prefix (alkyl/halogen), or the
    aromatic ring itself (phenyl or a heteroaromatic monocycle, via
    `name_branch`) when reached as one half of a
    `two_separate_rings_with_plain_aromatic_substituent` shape."""
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [
                name_branch(graph, root, atom, halogens, aromatic_atoms, mol=mol) for root in branch_roots
            ]
    return substituents


def _ring_name_from_substituents(ring_size, al_locant, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(ring_size, [], [], prefix, total_subs, "carbaldehyde", [al_locant])


def _ring_candidate_key(ring_size, al_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _ring_name_from_substituents(ring_size, al_locant, grouped)
    return al_locant, locant_set, citation_locants, name


def _name_ring_aldehyde(mol, ring_atoms, stereo=None, aromatic_atoms=frozenset()):
    """P-66.6.1.1.3: "The suffix 'carbaldehyde' is used when the -CHO
    group is attached to a carbon atom of a ring" -- e.g.
    'cyclohexanecarbaldehyde (PIN)'. Unlike the chain-parent 'al' suffix
    above (where the -CHO carbon is always the parent's own C1), here the
    ring itself is the parent hydride and the -CHO carbon is a
    substituent atom hanging directly off one ring carbon, mirroring
    `_carboxylic_acid.py`'s `_name_ring_carboxylic_acid` construction.

    A saturated monocyclic all-carbon ring with exactly one -CHO hanging
    directly off one ring atom (no intervening chain carbon, no
    standalone hydroxyl, no ring unsaturation), plus any number of
    substituents (alkyl/halogen) on any ring atom, including the -CHO-
    bearing atom itself.

    `stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- every stereocenter must lie on the ring
    itself (a stereocenter on the -CHO substituent branch is out of scope
    for now).

    `aromatic_atoms`: when the aldehyde-bearing ring reaches this function
    as one half of a `two_separate_rings_with_plain_aromatic_substituent`
    shape (mirroring `_nitrile.py`'s #631), the other, aromatic ring's
    atoms -- cited as a plain substituent (phenyl or a heteroaromatic
    monocycle) via `name_branch` the same way an ordinary alkyl ring
    substituent already is. Empty by default, so the original single-ring
    dispatch is unchanged."""
    aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol, aromatic_ring_atoms=aromatic_atoms)
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a ring aldehyde is not "
            "supported yet"
        )
    if len(aldehydes) != 1:
        raise UnsupportedStructure(
            "more than one aldehyde group alongside a ring is not "
            "supported yet"
        )
    (aldehyde_oxygen,) = aldehydes

    graph = adjacency(mol)
    (aldehyde_carbon,) = graph[aldehyde_oxygen]
    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] != aldehyde_oxygen
        and b[1] != aldehyde_oxygen
        and b[0] in ring_atoms
        and b[1] in ring_atoms
        and not (b[0] in aromatic_atoms and b[1] in aromatic_atoms)
    ]
    if all_non_single:
        raise UnsupportedStructure(
            "an unsaturated ring alongside an aldehyde substituent is not "
            "supported yet (see P-31.1.3)"
        )

    ring_neighbors = [n for n in graph[aldehyde_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "an aldehyde not directly attached to a single ring atom is "
            "not supported yet"
        )
    (ring_atom,) = ring_neighbors

    if stereo is not None and any(atom not in ring_atoms for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on the -CHO substituent branch rather than "
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
            al_locant = position_of[ring_atom]
            substituents = _ring_substituents(
                graph, candidate, halogens, {aldehyde_carbon}, mol=mol, aromatic_atoms=aromatic_atoms
            )
            key = _ring_candidate_key(ring_size, al_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _benzaldehyde_name_from_substituents(grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if total_subs == 0:
        # 'benzaldehyde' is a fully retained name (P-66.6.1.1.3) -- unlike
        # 'cyclohexanecarbaldehyde', there is no locant position to even
        # omit.
        return "benzaldehyde"
    return f"{format_substituent_prefixes(grouped)}benzaldehyde"


def _benzaldehyde_candidate_key(al_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzaldehyde_name_from_substituents(grouped)
    return al_locant, locant_set, citation_locants, name


def _name_benzaldehyde(mol, ring_atoms, exempt_atoms=None):
    """P-66.6.1.1.3: 'benzaldehyde' is a retained name that is itself the
    PIN for a -CHO hanging directly off one carbon of an otherwise-plain
    (or substituted) benzene ring -- e.g. 'benzaldehyde' (PubChem CID
    240), '2-methylbenzaldehyde' (CID 10722), '4-methylbenzaldehyde' (CID
    7725). Structurally identical to `_name_ring_aldehyde`'s saturated-
    ring case, mirroring `_carboxylic_acid.py`'s `_name_benzoic_acid`,
    with the retained name 'benzaldehyde' replacing 'cyclo' + alkane_name
    + 'carbaldehyde' as the whole suffix unit (no locant is ever cited
    for the -CHO position itself)."""
    aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol, aromatic_ring_atoms=exempt_atoms or ring_atoms)
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside benzaldehyde is not "
            "supported yet"
        )
    if len(aldehydes) != 1:
        raise UnsupportedStructure(
            "more than one aldehyde group alongside a benzene ring is not "
            "supported yet"
        )
    (aldehyde_oxygen,) = aldehydes

    graph = adjacency(mol)
    (aldehyde_carbon,) = graph[aldehyde_oxygen]
    ring_neighbors = [n for n in graph[aldehyde_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "an aldehyde not directly attached to a single ring atom is "
            "not supported yet"
        )
    (ring_atom,) = ring_neighbors
    other_ring_atom_branches = [n for n in graph[ring_atom] if n not in ring_atoms and n != aldehyde_carbon]
    if other_ring_atom_branches:
        raise UnsupportedStructure(
            "a substituent on the same ring atom as the aldehyde is not "
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
            al_locant = position_of[ring_atom]
            substituents = _ring_substituents(graph, candidate, halogens, {aldehyde_carbon}, mol=mol)
            key = _benzaldehyde_candidate_key(al_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    return best_name


def _name_phenyl_chain_aldehyde(mol, ring_atoms):
    """Name an aldehyde whose -CHO lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted benzene
    ring -- e.g. 3-phenylpropanal. The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring
    `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid`. Narrower
    than the acyclic path above: exactly one -CHO, no coexisting standalone
    hydroxyl, no chain unsaturation, and no specified stereocenter -- each
    is a separate follow-up (see
    tasks/phenyl-substituent-on-aldehyde-chain.md's scope note)."""
    aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol, aromatic_ring_atoms=ring_atoms)
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a benzene-ring-substituent "
            "aldehyde chain is not supported yet"
        )
    if len(aldehydes) != 1:
        raise UnsupportedStructure(
            "more than one aldehyde group alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "aldehyde chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in aldehydes and b[1] not in aldehydes and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "aldehyde chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_chain_attachments_with_halogens(graph, rings, set(), halogens)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain aldehyde is not "
            "supported yet"
        )
    (aldehyde_oxygen,) = aldehydes
    (aldehyde_carbon,) = graph[aldehyde_oxygen]
    chain, branches = longest_branched_chain(graph, aldehyde_carbon, ring_atoms, aldehydes, halogens=halogen_substituents(mol))

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, 1, [], [], grouped)


def name_aldehyde(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        union = set().union(*aromatic_rings)
        aldehydes, _ = _validate_and_collect_aldehydes(mol, aromatic_ring_atoms=union)
        anchors = [next(iter(adjacency(mol)[o])) for o in aldehydes]
        host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
        if host is not None:
            return _name_benzaldehyde(mol, host, union)
        return _name_phenyl_chain_aldehyde(mol, union)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        is_benzene = is_plain_benzene_ring(mol, ring_atoms)
        is_heteroaromatic = not is_benzene and (
            heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None
        )
        if is_benzene or is_heteroaromatic:
            aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol, aromatic_ring_atoms=ring_atoms)
            if not hydroxyls and len(aldehydes) == 1:
                (only_oxygen,) = aldehydes
                graph = adjacency(mol)
                (only_carbon,) = graph[only_oxygen]
                ring_neighbors = [n for n in graph[only_carbon] if n in ring_atoms]
                if len(ring_neighbors) == 1:
                    if is_benzene:
                        return _name_benzaldehyde(mol, ring_atoms)
                    raise UnsupportedStructure(
                        "an aldehyde directly attached to a heteroaromatic "
                        "ring is not supported yet (only a plain benzene "
                        "ring's 'benzaldehyde' direct-attachment path, or a "
                        "heteroaromatic ring cited as a chain substituent, "
                        "is)"
                    )
            return _name_phenyl_chain_aldehyde(mol, ring_atoms)
        aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol)
        if not hydroxyls and len(aldehydes) == 1:
            (only_oxygen,) = aldehydes
            graph = adjacency(mol)
            (only_carbon,) = graph[only_oxygen]
            ring_neighbors = [n for n in graph[only_carbon] if n in ring_atoms]
            if len(ring_neighbors) == 1:
                return _name_ring_aldehyde(mol, ring_atoms, specified_stereocenters(mol))

    # Two separate simple monocycles joined by one direct bond, one a plain
    # benzo/heteroaromatic ring with no substituent of its own (mirroring
    # `_nitrile.py`'s #631) -- e.g. 2-phenylcyclohexane-1-carbaldehyde --
    # reuses `_name_ring_aldehyde` with the aromatic ring's atoms passed
    # through as an exemption/substituent. Narrower than that: only the
    # "every aldehyde is ring-borne" split is handled below (matching this
    # shape's only verified real structures), mirroring #622/#624/#628/
    # #631's identical scoping decision.
    if ring_info.NumRings() == 2:
        aromatic_shape = two_separate_rings_with_plain_aromatic_substituent(mol, adjacency(mol))
        if aromatic_shape is not None:
            ring_atoms, aromatic_atoms, _, _ = aromatic_shape
            aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol, aromatic_ring_atoms=aromatic_atoms)
            if hydroxyls:
                raise UnsupportedStructure(
                    "a standalone hydroxyl alongside this two-ring "
                    "aromatic-substituent shape is not supported yet"
                )
            if len(aldehydes) != 1:
                raise UnsupportedStructure(
                    "more than one aldehyde group alongside this two-ring "
                    "aromatic-substituent shape is not supported yet"
                )
            (only_oxygen,) = aldehydes
            graph = adjacency(mol)
            (only_carbon,) = graph[only_oxygen]
            ring_neighbors = [n for n in graph[only_carbon] if n in ring_atoms]
            if len(ring_neighbors) != 1:
                raise UnsupportedStructure(
                    "an aldehyde on a chain hanging off the non-aromatic "
                    "ring, with the ring itself bearing no aldehyde of its "
                    "own, alongside this two-ring aromatic-substituent "
                    "shape is not supported yet"
                )
            return _name_ring_aldehyde(
                mol, ring_atoms, specified_stereocenters(mol), aromatic_atoms=aromatic_atoms
            )

    aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol)
    stereo = specified_stereo_elements(mol)
    graph = adjacency(mol)
    # Exclude each C=O carbonyl bond itself: `non_single_bonds` reports it as
    # order 2.0 same as a C=C, but it isn't a chain 'ene' bond (one endpoint
    # is the aldehyde oxygen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in aldehydes and b[1] not in aldehydes]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior carbonyl form and "
                "is out of scope for this module (P-31.1.4.2.4)"
            )

    if carbon_on_ring(mol, graph, [c for o in aldehydes for c in graph[o]]):
        raise UnsupportedStructure(
            "this ring shape alongside an aldehyde (more than one -CHO on "
            "the ring, other ring unsaturation, a standalone hydroxyl, or "
            "a ring other than a single saturated monocyclic/benzene one) "
            "is out of scope for this module's 'carbaldehyde' suffix path "
            "(P-33.3.1.2)"
        )
    return _name_acyclic_aldehyde(mol, aldehydes, hydroxyls, bonds, stereo)
