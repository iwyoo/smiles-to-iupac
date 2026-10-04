"""Naming of carboxylic acids (the '-oic acid' suffix, -COOH) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.1.1, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-65.1.1; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'oic acid' is the
  preselected suffix for -COOH, ranked senior to every other suffix this
  project handles (ester, amide, nitrile, aldehyde, ketone, alcohol, amine).
  Since it outranks all of them, a coexisting carbonyl (aldehyde/ketone-
  shaped) is still rejected as an unresolved competition; but a coexisting
  standalone hydroxyl (-OH), being junior to '-oic acid', is cited as the
  'hydroxy' substituent prefix instead (P-41), e.g. 'OC(=O)CCO' ->
  '3-hydroxypropanoic acid'. Any other heteroatom, or any oxygen not part of
  a full -COOH pattern or a standalone hydroxyl, is simply rejected.
- P-65.1.1.2: a -COOH substituent on a ring is named with the separate
  'carboxylic acid' suffix (ring name + 'carboxylic acid'), a different
  construction from the chain-terminal 'oic acid' suffix this module
  builds -- out of scope here (acyclic-chain-parent only). A ring may
  still appear as a plain 'phenyl' substituent on the chain instead (see
  `_name_phenyl_chain_carboxylic_acid` below) -- the narrowest slice of
  `tasks/aromatic-ring-substituent-parent-selection.md`'s largest single
  real-data coverage gap: a single, otherwise-unsubstituted benzene ring
  hanging off one end of an unbranched chain whose far end carries the
  sole -COOH, e.g. 'c1ccccc1CC(=O)O' -> '2-phenylethanoic acid' -- the
  systematic stem, not PubChem's retained 'phenylacetic acid' (PubChem CID
  999 for the structure match only; this module's own already-established
  convention already prefers the systematic stem once substituted, e.g.
  'ClCC(=O)O' -> '2-chloroethanoic acid', not '2-chloroacetic acid' --
  P-65.1.1 retains 'acetic acid' as PIN only for unsubstituted CH3COOH
  itself). A -COOH directly on a saturated monocyclic all-carbon ring is
  supported (see `_name_ring_carboxylic_acid`), including other ring
  substituents (alkyl/halogen), e.g. '4-methylcyclohexane-1-carboxylic
  acid' (PubChem CID 20330). A -COOH directly on a single benzene ring
  carbon is also supported (see `_name_benzoic_acid`): 'benzoic acid'
  (PIN, P-65.1.1's own retained name, CID 243) with or without other
  ring substituents, e.g. '2-methylbenzoic acid' (CID 8373); the
  retained name itself carries no locant for the -COOH position, unlike
  the cycloalkane case. Ring unsaturation outside a benzene ring, a
  standalone hydroxyl, a second -COOH, a polycyclic/naphthalene ring, or
  an intervening chain carbon between a saturated ring and the -COOH
  carbon remain out of scope and still raise `UnsupportedStructure`.
- A -COOH carbon is always a chain terminus: after its carbonyl (=O) and
  hydroxyl (-OH) oxygens, it has room for at most one more substituent,
  which must be another chain carbon (or nothing, for formic acid,
  H-COOH). Unlike `_alcohol.py`/`_ketone.py`, this means the parent chain's
  numbering is never a locant-minimization choice for the -COOH group
  itself: whichever chain end carries a -COOH carbon simply becomes C1
  (P-65.1.1, mirroring the aldehyde '-al' suffix's own terminal-carbon
  fixing).
- P-14.3.3: the suffix locant for -COOH (and 'dioic acid' for two -COOH
  groups) is never cited, since it is always fully determined by the -COOH
  carbon(s) being chain terminus/termini with no other possible position,
  e.g. 'propanoic acid' (not 'propan-1-oic acid'), 'hexanedioic acid' (not
  'hexane-1,6-dioic acid'). Confirmed against PubChem. Only 'ene'/'yne'
  locants (P-31.1.1.1-.2, same mechanics as `_alcohol.py`) are ever cited in
  this module's suffix, e.g. 'but-2-enoic acid', 'but-2-enedioic acid'
  (fumaric/maleic acid's PIN).
- P-44.1.1 (Chapter P-4): the senior parent structure captures the maximum
  number of principal characteristic groups (here, -COOH). As in
  `_alcohol.py`, this module only accepts a candidate chain that is
  simultaneously one of the longest carbon chains (P-44.3.2) and carries
  every -COOH-bearing carbon in the molecule; anything requiring a shorter
  principal chain to capture more -COOH groups is out of scope.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -oic acid suffix, reusing `halogen_substituents`/
  `format_substituent_prefixes` unchanged.
- P-91.3/P-92: a molecule
  with one or more *specified* tetrahedral stereocenters -- every one on
  the principal chain itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
  ascending locant order, e.g. '(2R)-2-chloropropanoic acid',
  '(2S,3S)-2-chloro-3-hydroxybutanoic acid' (a Blue Book worked example).
  CIP computation is delegated entirely to RDKit (`_common.specified_stereocenters`);
  this module only formats the result, same pattern as `_alcohol.py`. A
  stereocenter on a substituent branch, mixed with an unspecified one, or
  alongside E/Z double-bond stereo remains out of scope (raises
  `UnsupportedStructure`).

Explicitly out of scope (raise `UnsupportedStructure`):
- A ring bearing anything other than the single narrow shape
  `_name_ring_carboxylic_acid` supports (see above): unsaturated,
  aromatic, or polycyclic rings, a ring reached through an intervening
  chain carbon, or more than one -COOH group.
- Any oxygen that isn't part of a full -COOH pattern on some carbon, or a
  standalone hydroxyl on a carbon with no carbonyl (a lone carbonyl with no
  matching hydroxyl indicates an unresolved aldehyde/ketone competition; an
  ether or any other oxygen shape).
- A -COOH carbon with more than one carbon neighbor (impossible for a
  genuine carboxyl carbon, since its remaining two bonds are already the
  carbonyl and hydroxyl oxygens; guarded defensively).
- A standalone hydroxyl on a carbon that is also part of a C=C/C#C bond (an
  enol, tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
- Any other heteroatom (N, S, ...).
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
)
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_chain

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def has_carboxylic_acid_shape(mol) -> bool:
    """True if some carbon carries both a doubly-bonded, monovalent oxygen
    and a singly-bonded, monovalent, one-H oxygen (a -COOH pattern),
    regardless of whether the rest of the molecule is in scope. Used by
    `core.py` to route ahead of the ketone/alcohol dispatch, since a -COOH
    carbon would otherwise look aldehyde- or alcohol-shaped to those
    modules."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        has_carbonyl = any(
            o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            for o in oxygens
        )
        has_hydroxyl = any(
            o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and o.GetTotalNumHs() == 1
            for o in oxygens
        )
        if has_carbonyl and has_hydroxyl:
            return True
    return False


def _validate_and_collect_carboxyls(mol, aromatic_ring_atoms=frozenset(), extra_accounted_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (carboxyl_carbons, carboxyl_oxygens, extra_hydroxyls): the set
    of -COOH carbon atom indices, the set of all their carbonyl/hydroxyl
    oxygen atom indices (excluded from substituent detection), and the set of
    any coexisting standalone alcohol hydroxyl-oxygen atom indices. A
    standalone hydroxyl is junior to 'oic acid' in Table 3.3's suffix
    seniority order, so it is cited as the 'hydroxy' substituent prefix
    instead of competing for the suffix (P-41).

    `aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene or heteroaromatic-monocycle ring (P-29.3.4.1) with exactly one
    exocyclic attachment -- exempted from both the atomic-number allowlist
    and the aromatic-atom rejection below (including the ring's own
    heteroatom, if it has one) so `name_carboxylic_acid`'s ring-substituent
    path (see `_name_phenyl_chain_carboxylic_acid`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged.

    `extra_accounted_atoms`: atom indices a coexisting-group pairwise
    module (e.g. `_carboxylic_acid_sulfinic_acid.py`) has already
    independently validated as belonging to its own demoted junior group
    (its heteroatom(s) plus their own oxygens) -- skipped entirely here
    (not subject to the plain-carboxylic-acid atom-type allowlist below),
    reused via `_coexisting_groups.name_via_senior_phenyl_chain`. Empty by
    default, so every other caller's behavior is unchanged."""
    has_carbon = False
    carbonyls_by_carbon = {}
    hydroxyls_by_carbon = {}
    for atom in mol.GetAtoms():
        if atom.GetIdx() in extra_accounted_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS and atom.GetIdx() not in aromatic_ring_atoms:
            raise UnsupportedStructure(
                "heteroatoms other than carboxylic-acid oxygens (P-65.1.1) "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetIdx() in aromatic_ring_atoms and atomic_num != 6:
            continue
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether or ester) is out of scope; only oxygens forming "
                    "an isolated -COOH group are supported (P-65.1.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a -COOH oxygen must be attached to a carbon atom")
            bond_order = bond.GetBondTypeAsDouble()
            if bond_order == 2.0:
                carbonyls_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
            elif bond_order == 1.0 and atom.GetTotalNumHs() == 1:
                hydroxyls_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
            else:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "half of a -COOH group is out of scope for this module"
                )
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

    carboxyl_carbons = set()
    carboxyl_oxygens = set()
    extra_hydroxyls = set()
    for carbon_idx in set(carbonyls_by_carbon) | set(hydroxyls_by_carbon):
        carbonyls = carbonyls_by_carbon.get(carbon_idx, [])
        hydroxyls = hydroxyls_by_carbon.get(carbon_idx, [])
        if not carbonyls and len(hydroxyls) == 1:
            # A hydroxyl-only carbon is a standalone alcohol (-OH), junior to
            # -COOH in Table 3.3 (see docstring) - cited as a 'hydroxy'
            # prefix rather than rejected outright.
            extra_hydroxyls.add(hydroxyls[0])
            continue
        if len(carbonyls) != 1 or len(hydroxyls) != 1:
            raise UnsupportedStructure(
                "an oxygen pattern that isn't exactly one carbonyl and one "
                "hydroxyl oxygen on the same carbon, or a standalone "
                "hydroxyl, is a more/less senior characteristic group than a "
                "plain carboxylic acid (Table 3.3), which this module does "
                "not attempt to disambiguate"
            )
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "a -COOH carbon with more than one carbon neighbor is not a "
                "valid carboxyl carbon"
            )
        carboxyl_carbons.add(carbon_idx)
        carboxyl_oxygens.add(carbonyls[0])
        carboxyl_oxygens.add(hydroxyls[0])

    if not carboxyl_carbons:
        raise UnsupportedStructure(
            "no carboxylic acid (-COOH) group found; this module only "
            "handles carboxylic acids"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return carboxyl_carbons, carboxyl_oxygens, extra_hydroxyls


def _name_from_substituents(chain_length, acid_count, ene_locants, yne_locants, grouped):
    # P-14.3.4.2(a): a mononuclear parent's substituent locant is always
    # '1' and never cited.
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    return (
        prefix
        + name_from_substituents(chain_length, ene_locants, yne_locants, multiplied_word(acid_count, "oic"))
        + " acid"
    )


def _candidate_key(chain_length, acid_count, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.10 (ene/yne locants) ahead of P-45.2
    (substituent-prefix locants), most-preferred first. The -COOH group's
    own locant isn't part of this key: candidates are pre-filtered so a
    -COOH carbon always sits at C1 (see `_name_acyclic_carboxylic_acid`)."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, acid_count, ene_locants, yne_locants, grouped)
    return (
        (
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _best_acyclic_carboxylic_acid_candidate(
    mol,
    carboxyl_carbons,
    carboxyl_oxygens,
    hydroxyls,
    bonds,
    stereo=None,
    extra_names=None,
    required_atoms=frozenset(),
):
    """(best_name, best_position_of) -- the winning chain numbering and its
    fully-formatted name, factored out of `_name_acyclic_carboxylic_acid`
    (which just adds the stereo-descriptor prefix on top) so
    `_isotope_carboxylic_acid.py` can reuse the identical numbering
    decision, mirroring `_alcohol.py`'s/`_ketone.py`'s identical factoring.

    `extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group demoted to a substituent prefix by
    `_seniority.senior_class` (e.g. a demoted amine's 'amino'), reused by
    `_coexisting_groups.py` so a pairwise module doesn't have to
    reimplement this function's chain search/candidate selection. `None`
    keeps the original hydroxyl-only behavior unchanged. `required_atoms`:
    additional carbon atoms (e.g. every demoted amine's own carbon
    neighbor) that a candidate chain must also carry -- empty by default
    so existing callers are unaffected."""
    graph = adjacency(mol)
    halogens = {
        **halogen_substituents(mol),
        **{o: "hydroxy" for o in hydroxyls},
        **(extra_names or {}),
    }
    chains = all_chains(carbon_adjacency(mol))
    acid_count = len(carboxyl_carbons)
    stereo_atoms = stereo_element_atoms(mol, stereo) if stereo is not None else []

    eligible = []
    for chain in chains:
        chain_set = set(chain)
        if not carboxyl_carbons <= chain_set:
            continue
        if not required_atoms <= chain_set:
            continue
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            carboxyl_carbons <= set(c)
            and required_atoms <= set(c)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereo element on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every -COOH-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more -COOH groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] not in carboxyl_carbons:
                # A -COOH carbon must sit at C1 (see module docstring); a
                # direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, carboxyl_oxygens, mol=mol)
            key, name = _candidate_key(chain_length, acid_count, ene_locants, yne_locants, substituents)
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            key = (key, stereo_locant_rank(mol, stereo, position_of))
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    return best_name, best_position_of


def _name_acyclic_carboxylic_acid(
    mol,
    carboxyl_carbons,
    carboxyl_oxygens,
    hydroxyls,
    bonds,
    stereo=None,
    extra_names=None,
    required_atoms=frozenset(),
):
    """`stereo`: None, or a list of ("atom"/"bond", idx, "R"/"S"/"E"/"Z")
    from `specified_stereo_elements`/`specified_stereocenters` -- if given,
    only chain candidates that include every tetrahedral stereocenter are
    eligible (P-92: a stereocenter on a substituent branch rather than the
    principal chain is out of scope, mirroring `_alcohol.py`'s
    `_name_acyclic_alcohol`; a double-bond E/Z element's atoms are already
    required to lie on the chain via `bonds`, so no separate check is
    needed for those), and the winning candidate's own locants are used to
    format a "(<locant><R/S/E/Z>,...)-" prefix onto the name, ascending
    locant order (P-91.3, including when both kinds coexist) -- same
    mechanism as `_alcohol.py`, since a -COOH carbon's own fixed C1
    position (see module docstring) already decides numbering before
    stereo is even considered. See `_best_acyclic_carboxylic_acid_candidate`
    for `extra_names`/`required_atoms`."""
    best_name, best_position_of = _best_acyclic_carboxylic_acid_candidate(
        mol, carboxyl_carbons, carboxyl_oxygens, hydroxyls, bonds, stereo, extra_names, required_atoms
    )
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


def _name_phenyl_chain_carboxylic_acid(
    mol, ring_atoms, extra_names=None, required_atoms=frozenset(), extra_accounted_atoms=frozenset()
):
    """Name a carboxylic acid whose -COOH lies on a chain hanging off one
    atom of a benzene or plain heteroaromatic-monocycle ring (pyridine/
    furan/thiophene/pyrrole, P-29.3.4.1) -- e.g. 2-phenylethanoic acid, or
    2-(pyridin-3-yl)acetic acid (PubChem CID 108). The ring is cited as a
    'phenyl'/'pyridin-3-yl'/etc. (or, for a plain benzene ring whose other
    atoms each carry a single halogen/alkyl substituent, e.g.
    '4-chlorophenyl') substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_alcohol.py`'s `_name_ring_substituent_chain_alcohol` for a plain
    saturated ring. The chain itself may branch (`longest_branched_chain`
    picks one of the longest chains containing the -COOH carbon, per
    P-44.3.2, absorbing a branch into the parent chain whenever that makes
    it longer) -- e.g. ibuprofen's alpha-methyl becomes part of 'propanoic
    acid' rather than a substituent on a shorter 'ethanoic acid'
    (`CC(C)Cc1ccc(cc1)C(C)C(=O)O` -> '2-[4-(2-methylpropyl)phenyl]
    propanoic acid', PubChem CID 3672). Narrower than the acyclic path
    above: exactly one -COOH, no coexisting standalone hydroxyl, no
    chain unsaturation, no specified stereocenter, and (for a plain
    benzene ring) no non-halogen/non-alkyl ring substituent alongside the
    chain (a substituent on the heteroaromatic ring's own other atoms is
    out of scope for now, same restriction) -- each is a separate
    follow-up rather than being combined with the ring case in this first
    slice.

    `extra_names`/`required_atoms`/`extra_accounted_atoms`: same injection
    point as `_name_acyclic_carboxylic_acid`
    (`extra_accounted_atoms` mirrors that function's caller-side atom-type
    validation, since this function -- unlike the acyclic one -- still
    does its own), reused by `_coexisting_groups.
    name_via_senior_phenyl_chain` so a pairwise module doesn't have to
    reimplement this function's chain search/numbering logic --
    `None`/empty keeps the original behavior unchanged."""
    carboxyl_carbons, carboxyl_oxygens, hydroxyls = _validate_and_collect_carboxyls(
        mol, aromatic_ring_atoms=ring_atoms, extra_accounted_atoms=extra_accounted_atoms
    )
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a benzene-ring-substituent "
            "carboxylic acid chain is not supported yet"
        )
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "more than one carboxylic acid group alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "carboxylic acid chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in carboxyl_oxygens
        and b[1] not in carboxyl_oxygens
        and b[0] not in ring_atoms
        and b[1] not in ring_atoms
        and b[0] not in extra_accounted_atoms
        and b[1] not in extra_accounted_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "carboxylic acid chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **(extra_names or {})}
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_branch_attachments(mol, graph, rings, known=extra_names or ())
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain carboxylic acid is "
            "not supported yet"
        )
    (carboxyl_carbon,) = carboxyl_carbons
    chain, branches = longest_branched_chain(
        graph,
        carboxyl_carbon,
        ring_atoms,
        carboxyl_oxygens,
        halogens={**halogen_substituents(mol), **(extra_names or {})},
    )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a -COOH group directly attached to the benzene ring (no "
            "intervening chain carbon) uses the separate 'carboxylic "
            "acid' suffix construction (P-65.1.1.2), out of scope for "
            "this acyclic-chain-parent module"
        )
    if not required_atoms <= set(chain):
        raise UnsupportedStructure(
            "a required substituent atom outside the single unbranched "
            "chain hanging off the benzene ring is not supported yet"
        )

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, 1, [], [], grouped)


def _ring_substituents(graph, ring_order, halogens, excluded, mol=None):
    """{ring position -> [substituent name, ...]}, mirroring
    `_sulfonic_acid.py`'s identically-named helper -- every branch hanging
    off a ring atom other than the carboxylic acid carbon itself (in
    `excluded`) is a plain substituent prefix (alkyl/halogen)."""
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, carboxyl_locant, grouped, suffix="carboxylic acid"):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(ring_size, [], [], prefix, total_subs, suffix, [carboxyl_locant])


def _ring_candidate_key(ring_size, carboxyl_locant, substituents, suffix="carboxylic acid"):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _ring_name_from_substituents(ring_size, carboxyl_locant, grouped, suffix)
    return carboxyl_locant, locant_set, citation_locants, name


def _name_ring_attached_carboxyl(graph, ring_atoms, carboxyl_carbon, halogens, stereo=None, suffix="carboxylic acid", mol=None):
    """Shared ring-numbering-search kernel behind `_name_ring_carboxylic_
    acid` -- also reused by `_ester.py` (with `suffix="carboxylate"`) for
    the analogous ester-acyl case, since the numbering/substituent-citation
    logic is identical, only the trailing suffix word differs. The caller
    must have already validated: exactly one carboxyl-shaped carbon,
    directly attached to exactly one ring atom, no other substituent
    sharing that ring atom, and (if `stereo` is given) every stereocenter
    on the ring itself."""
    (ring_atom,) = [n for n in graph[carboxyl_carbon] if n in ring_atoms]
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            carboxyl_locant = position_of[ring_atom]
            substituents = _ring_substituents(graph, candidate, halogens, {carboxyl_carbon}, mol=mol)
            key = _ring_candidate_key(ring_size, carboxyl_locant, substituents, suffix)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_ring_carboxylic_acid(mol, ring_atoms, stereo=None):
    """P-65.1.2.2.2: "Carboxy groups attached to cyclic parent hydrides ...
    are always named by using the suffix 'carboxylic acid'" -- e.g.
    'cyclopentanecarboxylic acid (PIN)'. Unlike the chain-parent 'oic acid'
    suffix above (where the -COOH carbon is always the parent's own C1),
    here the ring itself is the parent hydride and the -COOH carbon is a
    substituent atom hanging directly off one ring carbon, mirroring
    `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid` construction but for
    an *exocyclic* suffix carbon rather than a ring-atom-attached one.

    A saturated monocyclic all-carbon ring with exactly one -COOH hanging
    directly off one ring atom (no intervening chain carbon, no standalone
    hydroxyl, no ring unsaturation), plus any number of substituents
    (alkyl/halogen) on any ring atom, including the -COOH-bearing atom
    itself -- e.g. 'cyclohexanecarboxylic acid' (the sole substituent's
    ring locant is P-14.3.3-omitted), '4-methylcyclohexane-1-carboxylic
    acid' (PubChem CID 20330), and '1-methylcyclohexane-1-carboxylic
    acid' (quaternary -COOH-bearing atom), mirroring `_sulfonic_acid.py`'s
    ring-numbering search.

    `stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- every stereocenter must lie on the ring
    itself (a stereocenter on the -COOH substituent branch is out of
    scope for now, mirroring `_sulfonic_acid.py`'s
    `_name_cyclic_sulfonic_acid` before its own branch-stereo extension)."""
    carboxyl_carbons, carboxyl_oxygens, extra_hydroxyls = _validate_and_collect_carboxyls(mol)
    if extra_hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a ring carboxylic acid is not "
            "supported yet"
        )
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "more than one carboxylic acid group alongside a ring is not "
            "supported yet"
        )
    (carboxyl_carbon,) = carboxyl_carbons

    all_non_single = non_single_bonds(mol)
    if any(
        a not in carboxyl_oxygens and b not in carboxyl_oxygens and a in ring_atoms and b in ring_atoms
        for a, b, _ in all_non_single
    ):
        raise UnsupportedStructure(
            "an unsaturated ring alongside a carboxylic acid substituent "
            "is not supported yet (see P-31.1.3)"
        )

    graph = adjacency(mol)
    ring_neighbors = [n for n in graph[carboxyl_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "a carboxylic acid not directly attached to a single ring atom "
            "is not supported yet"
        )
    (ring_atom,) = ring_neighbors

    if stereo is not None and any(atom not in ring_atoms for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on the -COOH substituent branch rather than "
            "the ring itself is not supported yet (see P-92)"
        )

    halogens = halogen_substituents(mol)
    return _name_ring_attached_carboxyl(graph, ring_atoms, carboxyl_carbon, halogens, stereo, mol=mol)


def _benzoic_acid_name_from_substituents(grouped, word="benzoic acid"):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if total_subs == 0:
        # 'benzoic acid' is a fully retained name (P-65.1.1) -- unlike
        # 'cyclohexanecarboxylic acid', there is no locant position to
        # even omit.
        return word
    return f"{format_substituent_prefixes(grouped)}{word}"


def _benzoic_acid_candidate_key(carboxyl_locant, substituents, word="benzoic acid"):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzoic_acid_name_from_substituents(grouped, word)
    return carboxyl_locant, locant_set, citation_locants, name


def _name_benzo_attached_carboxyl(graph, ring_atoms, carboxyl_carbon, halogens, word="benzoic acid", mol=None):
    """Shared ring-numbering-search kernel behind `_name_benzoic_acid` --
    also reused by `_ester.py` (with `word="benzoate"`) for the analogous
    ester-acyl case (P-65.6.3.2's retained-name-as-acyl-part axis). The
    caller must have already validated: exactly one carboxyl-shaped
    carbon, directly attached to exactly one ring atom of an otherwise-
    plain benzene ring."""
    (ring_atom,) = [n for n in graph[carboxyl_carbon] if n in ring_atoms]
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            carboxyl_locant = position_of[ring_atom]
            substituents = _ring_substituents(graph, candidate, halogens, {carboxyl_carbon}, mol=mol)
            key = _benzoic_acid_candidate_key(carboxyl_locant, substituents, word)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_benzoic_acid(mol, ring_atoms, exempt_atoms=None):
    """P-65.1.1: 'benzoic acid' is a retained name that is itself the PIN
    for a -COOH hanging directly off one carbon of an otherwise-plain (or
    substituted) benzene ring -- e.g. 'benzoic acid' (PubChem CID 243),
    '2-methylbenzoic acid' (CID 8373), '4-methylbenzoic acid' (CID 7470).
    Structurally identical to `_name_ring_carboxylic_acid`'s saturated-
    ring case (the -COOH carbon is an exocyclic substituent atom, not a
    ring atom itself) -- mirrors that function's ring-numbering search,
    with the retained name 'benzoic acid' replacing 'cyclo' + alkane_name
    + 'carboxylic acid' as the whole suffix unit (no locant is ever cited
    for the -COOH position itself, unlike the cycloalkane case, since
    'benzoic acid' carries no positional stem at all)."""
    carboxyl_carbons, carboxyl_oxygens, extra_hydroxyls = _validate_and_collect_carboxyls(
        mol, aromatic_ring_atoms=exempt_atoms or ring_atoms
    )
    if extra_hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside benzoic acid is not "
            "supported yet"
        )
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "more than one carboxylic acid group alongside a benzene ring "
            "is not supported yet"
        )
    (carboxyl_carbon,) = carboxyl_carbons

    graph = adjacency(mol)
    ring_neighbors = [n for n in graph[carboxyl_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "a carboxylic acid not directly attached to a single ring atom "
            "is not supported yet"
        )
    (ring_atom,) = ring_neighbors

    halogens = halogen_substituents(mol)
    return _name_benzo_attached_carboxyl(graph, ring_atoms, carboxyl_carbon, halogens, mol=mol)


def _has_carboxyl_directly_on_ring(mol, ring_atoms):
    """Cheap routing check (not full validation): True iff some ring atom
    has an exocyclic carbon neighbor shaped like a -COOH carbon (two
    oxygen neighbors) -- used only to decide between the benzoic-acid
    ring-parent path and the phenyl-chain path below."""
    for idx in ring_atoms:
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            if nbr.GetIdx() in ring_atoms or nbr.GetAtomicNum() != 6:
                continue
            oxygens = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 8]
            if len(oxygens) == 2:
                return True
    return False


def name_carboxylic_acid(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        union = set().union(*aromatic_rings)
        carboxyl_carbons, *_ = _validate_and_collect_carboxyls(mol, aromatic_ring_atoms=union)
        anchors = list(carboxyl_carbons)
        host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
        if host is not None:
            return _name_benzoic_acid(mol, host, union)
        return _name_phenyl_chain_carboxylic_acid(mol, union)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            if _has_carboxyl_directly_on_ring(mol, ring_atoms):
                return _name_benzoic_acid(mol, ring_atoms)
            return _name_phenyl_chain_carboxylic_acid(mol, ring_atoms)
        if heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None:
            if _has_carboxyl_directly_on_ring(mol, ring_atoms):
                raise UnsupportedStructure(
                    "a -COOH group directly on a heteroaromatic ring uses "
                    "a different suffix-naming construction, out of scope "
                    "for this module's chain-substituent path"
                )
            return _name_phenyl_chain_carboxylic_acid(mol, ring_atoms)
        if not any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
            if _has_carboxyl_directly_on_ring(mol, ring_atoms):
                return _name_ring_carboxylic_acid(mol, ring_atoms, specified_stereocenters(mol))
        else:
            raise UnsupportedStructure(
                "a -COOH group on/in a ring uses the separate 'carboxylic acid' "
                "suffix construction (P-65.1.2.2.2), out of scope for this "
                "acyclic-only module"
            )
    elif ring_info.NumRings() > 0:
        raise UnsupportedStructure(
            "a -COOH group on/in a ring uses the separate 'carboxylic acid' "
            "suffix construction (P-65.1.2.2.2), out of scope for this "
            "acyclic-only module"
        )
    carboxyl_carbons, carboxyl_oxygens, hydroxyls = _validate_and_collect_carboxyls(mol)
    stereo = specified_stereo_elements(mol)

    all_non_single = [
        b for b in non_single_bonds(mol) if b[0] not in carboxyl_oxygens and b[1] not in carboxyl_oxygens
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    graph = adjacency(mol)
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior carbonyl form and "
                "is out of scope for this module (P-31.1.4.2.4)"
            )

    return _name_acyclic_carboxylic_acid(mol, carboxyl_carbons, carboxyl_oxygens, hydroxyls, bonds, stereo)
