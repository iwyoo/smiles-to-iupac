"""Naming of primary amides (the '-amide' suffix, terminal -C(=O)NH2) on
acyclic saturated or unsaturated carbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-66.1, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-66.1; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'amide' is the preselected
  suffix for -CONH2, ranked junior to '-oic acid'/ester but senior to
  nitrile/aldehyde/ketone/alcohol/amine. A coexisting standalone hydroxyl
  (-OH), being junior to 'amide', is cited as the 'hydroxy' substituent
  prefix instead (P-41), mirroring `_aldehyde.py`/`_carboxylic_acid.py`. Any
  other coexisting carbonyl (aldehyde/ketone/-COOH-shaped) is rejected as an
  unresolved suffix-vs-suffix seniority competition, same as those modules.
- Like a -COOH carbon (`_carboxylic_acid.py`), an amide carbon is always a
  chain terminus: after its carbonyl (=O) and amide nitrogen, it has room
  for at most one more substituent, which must be another chain carbon (or
  nothing, for formamide, HCONH2). So the parent chain's
  numbering is never a locant-minimization choice for the amide group
  itself: whichever chain end carries the amide carbon simply becomes C1.
- P-14.3.3: the suffix locant for -CONH2 is never cited, since it is always
  fully determined by the amide carbon being a chain terminus with no other
  possible position, e.g. 'propanamide' (not 'propan-1-amide').
- P-31.0/P-31.1.1.1-.2: construction of the 'ene'/'yne' portion of a combined
  unsaturated-amide name reuses the same mechanics as `_aldehyde.py`/
  `_carboxylic_acid.py` (ending replaces 'ane' entirely, multiplying
  prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before 'yne', euphonic
  stem 'a' before a multiplied ending); 'amide' is appended directly onto
  the last such segment (eliding a trailing vowel the same way 'ene'+'al'
  does), since it never carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'amide' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.
- The amide nitrogen may carry zero, one, or two plain,
  unsubstituted, saturated, acyclic alkyl substituents (branched or
  unbranched), each cited as its own 'N-'-prefixed substituent directly
  ahead of the acyl stem, in alphanumerical order (P-14.5.2, ignoring
  italicized prefixes like 'tert-' via `alpha_sort_key`), with a 'di'
  multiplying prefix (and a single shared 'N,N-' pair) when both are
  identical -- exactly `_carbamate.py`'s/`_urea.py`'s own N,N-
  disubstitution extension. Each N-substituent's own name is built with
  `name_branch` (P-29 PIN style, fixed project-wide by PR #237; mirrors
  `_carbamate.py`'s/`_urea.py`'s identical fix, PR #328/#332), e.g.
  'N-tert-butylacetamide' (CID 12985, a retained non-compound name). A
  *compound* N-substituent (has its own locant) is always parenthesized
  -- 'N-(propan-2-yl)acetamide', not PubChem's own raw
  'N-propan-2-ylacetamide' (CID 136874), same correction as `_urea.py`
  (see that module's docstring for the Blue Book citations). A multiplied
  identical-pair name is likewise parenthesized only when compound, e.g.
  'N,N-di(propan-2-yl)acetamide' (CID 69797) vs. 'N,N-ditert-butylacetamide'
  (CID 18999412). Confirmed via PubChem: 'N-methylacetamide' (CC(=O)NC),
  'N,N-dimethylacetamide' (CC(=O)N(C)C), 'N-ethyl-N-methylacetamide'
  (CC(=O)N(C)CC), 'N-tert-butyl-N-ethylacetamide' (CID 54197906,
  alphabetized ignoring 'tert-').

- P-91.3/P-92: a molecule with one
  or more *specified* tetrahedral stereocenters -- every one on the
  principal chain itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-"
  prefix, ascending locant order, cited outermost (ahead of any N-alkyl
  prefix, since a stereodescriptor always sits at the very front of the
  complete name), e.g. '(2R)-2-methylbutanamide',
  '(2R)-N-ethyl-2-methylbutanamide'. The amide nitrogen itself (planar,
  sp2) is never a stereocenter, so this is unconditional like
  `_carboxylic_acid.py`/`_aldehyde.py`. While wiring this in, a
  pre-existing gap surfaced: an N-substituent bearing its own hydroxyl
  (invisible to the carbon-only chain walk that measures N-substituent
  length) used to be silently accepted and misnamed as if it were plain
  alkyl -- now explicitly rejected too.

Explicitly out of scope (raise `UnsupportedStructure`), per the task's
first-pass scope:
- An N-substituent that is unsaturated or ring-bearing (a branched but
  otherwise plain saturated acyclic N-substituent is supported, see
  above), except for one narrow case: a plain, unsubstituted benzene ring
  hanging directly off the amide nitrogen itself (an anilide, e.g.
  'N-phenylacetamide', PubChem confirms `CC(=O)Nc1ccccc1` ->
  'N-phenylacetamide') -- mirrors the identical 'phenyl N-substituent'
  pattern already wired into `_urea.py`/`_thiourea.py`/`_guanidine.py`/
  `_selenourea.py`/`_tellurourea.py`/`_carbamate.py`. Only a lone phenyl
  (no second substituent sharing that nitrogen, no substituted phenyl
  like 4-hydroxyphenyl) is supported; each is a separate follow-up.
- A true lactam (the carbonyl carbon itself is a ring atom) - handled by
  `_ketone.py`'s hetero-ring ketone path, not this module. An
  N-unsubstituted -CONH2 hanging as an *exocyclic* substituent directly
  off one ring carbon (P-66.1.1.1.1.3's 'carboxamide' suffix) is handled
  by `_name_ring_amide`/`_name_benzamide` for a single saturated
  monocyclic or benzene ring only - N-substitution, multiple amides, a
  standalone hydroxyl, or ring unsaturation alongside it are still out
  of scope.
- More than one amide group in the same molecule (a diamide) - deferred
  entirely, along with any other multiple-principal-characteristic-group
  combination.
- Any oxygen that isn't the amide's own carbonyl oxygen or a standalone
  hydroxyl (ethers, esters, -COOH, and any oxygen bonded to more than one
  heavy atom).
- A carbonyl carbon with other than exactly one nitrogen neighbor, or more
  than one carbon neighbor (a ketone-shaped carbon is not an amide carbon).
- An aromatic carbonyl carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory, except for one narrow case: a
  primary amide's chain hanging off a single plain, unsubstituted benzene
  ring with no other substituent on the ring (`_name_phenyl_chain_amide`,
  e.g. '3-phenylpropanamide'), mirroring `_aldehyde.py`/
  `_carboxylic_acid.py`'s identical benzene-ring-substituent path. Narrower
  than the acyclic path: no N-alkyl substitution, no coexisting standalone
  hydroxyl, no chain unsaturation, no specified stereocenter, and no
  substituted benzene/naphthalene - each a separate follow-up.
- Any other heteroatom (S, ...).
- A hydroxyl on a carbon that is also part of a C=C/C#C bond (an enol,
  tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    all_chains,
    bfs,
    carbon_adjacency,
    carbon_on_ring,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    lowest_locant_set,
    most_multiple_bonds,
    name_from_substituents,
    non_single_bonds,
    ring_chain_attachment,
    ring_branch_attachments,
    ring_hosting_anchors,
    separate_aromatic_monocycles,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._substituents import (
    format_substituent_prefixes,
    name_branch,
    substituents_for_chain,
)

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _is_carbonyl_carbon(mol, carbon_atom):
    return any(
        o.GetAtomicNum() == 8
        and o.GetDegree() == 1
        and mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        for o in carbon_atom.GetNeighbors()
    )


def _is_plain_hydroxyl_oxygen(mol, oxygen_atom, other_atom_idx):
    """True if `oxygen_atom` is a bare -OH: degree 1, singly bonded to
    `other_atom_idx`, and exactly one implicit/explicit hydrogen. Shared by
    `has_amide_shape` and `_validate_and_collect_amide` to recognize a
    hydroxamic acid's N-hydroxy oxygen (P-66.1.1.3.2) the same way as any
    other plain standalone hydroxyl in this module."""
    return (
        oxygen_atom.GetAtomicNum() == 8
        and oxygen_atom.GetDegree() == 1
        and oxygen_atom.GetTotalNumHs() == 1
        and mol.GetBondBetweenAtoms(other_atom_idx, oxygen_atom.GetIdx()).GetBondTypeAsDouble() == 1.0
    )


def has_amide_shape(mol) -> bool:
    """True if some carbon carries a doubly-bonded, monovalent carbonyl
    oxygen and a singly-bonded nitrogen with 0-2 carbon substituents (or a
    single plain hydroxyl substituent, P-66.1.1.3.2's N-hydroxy amide/
    hydroxamic acid shape), no other heavy-atom neighbor, and no *other*
    carbonyl-carbon neighbor (a -CON(R)(R') pattern, R/R' either H, an
    unbranched alkyl carbon, or -OH), regardless of whether the rest of the
    molecule is in scope. Used by `core.py` to route ahead of the aldehyde/
    ketone dispatch, since an amide carbon would otherwise look aldehyde-
    shaped to those modules (both have exactly one carbon neighbor besides
    the carbonyl, and that count ignores a nitrogen substituent -- the same
    routing trap `core.py`'s own `_is_aldehyde_shaped` comment already names
    directly). A nitrogen bonded to two carbonyl carbons (a symmetric imide)
    is excluded here so `core.py`'s later `has_imide_shape` check still gets
    a chance at it."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if not _is_carbonyl_carbon(mol, atom):
            continue
        has_amide_n = any(
            n.GetAtomicNum() == 7
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and all(
                nn.GetAtomicNum() == 6 or _is_plain_hydroxyl_oxygen(mol, nn, n.GetIdx())
                for nn in n.GetNeighbors()
                if nn.GetIdx() != atom.GetIdx()
            )
            and sum(1 for nn in n.GetNeighbors() if _is_carbonyl_carbon(mol, nn)) == 1
            for n in atom.GetNeighbors()
        )
        if has_amide_n:
            return True
    return False


def _validate_and_collect_amide(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons,
    hydroxyls, n_hydroxy_oxygen): the single -CON(R)(R') carbon/oxygen/
    nitrogen atom indices, a tuple of 0-2 N-alkyl substituent carbon
    indices, the set of any coexisting standalone (chain) hydroxyl-oxygen
    atom indices, and the amide nitrogen's own hydroxyl-oxygen atom index
    if it has one (a hydroxamic acid, P-66.1.1.3.2) else None.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_amide`'s benzene-ring-substituent
    path (see `_name_phenyl_chain_amide`) can reuse this same validation for
    the rest of the molecule. Empty by default, so every other caller's
    behavior is unchanged."""
    has_carbon = False
    amide_carbons = set()
    amide_oxygen_by_carbon = {}
    amide_nitrogen_by_carbon = {}
    n_alkyl_carbons_by_nitrogen = {}
    n_hydroxy_by_nitrogen = {}
    hydroxyls = set()
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a primary amide's oxygen/nitrogen "
                "(P-66.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
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
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether or ester) is out of scope; only an isolated "
                    "amide carbonyl or a standalone hydroxyl is supported "
                    "(P-66.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() == 7:
                # A hydroxamic acid's N-hydroxy oxygen (P-66.1.1.3.2) --
                # the nitrogen branch below validates and collects it, not
                # this generic chain/ring-hydroxyl path.
                continue
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("an amide/hydroxyl oxygen must be attached to a carbon atom")
            bond_order = bond.GetBondTypeAsDouble()
            if bond_order == 1.0 and atom.GetTotalNumHs() == 1:
                hydroxyls.add(atom.GetIdx())
                continue
            if bond_order != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring is out of scope for "
                    "this module"
                )
            amide_oxygen_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
        elif atomic_num == 7:
            neighbors = list(atom.GetNeighbors())
            non_carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() != 6]
            if len(non_carbon_neighbors) > 1 or any(
                not _is_plain_hydroxyl_oxygen(mol, n, atom.GetIdx()) for n in non_carbon_neighbors
            ):
                raise UnsupportedStructure(
                    "an amide nitrogen bonded to anything other than "
                    "carbon, or a single plain hydroxyl (P-66.1.1.3.2), is "
                    "out of scope for this module"
                )
            n_hydroxy_oxygen = non_carbon_neighbors[0].GetIdx() if non_carbon_neighbors else None
            if len(neighbors) > 3:
                raise UnsupportedStructure(
                    "an amide nitrogen with more than two substituents "
                    "besides its carbonyl carbon is not a valid amide "
                    "nitrogen"
                )
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() != 1.0:
                    raise UnsupportedStructure(
                        "an amide nitrogen must be singly bonded to all its neighbors"
                    )

            carbonyl_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6 and _is_carbonyl_carbon(mol, n)]
            if len(carbonyl_neighbors) != 1:
                raise UnsupportedStructure(
                    "a nitrogen bonded to zero or multiple carbonyl carbons "
                    "is not a valid amide nitrogen for this module"
                )
            (carbonyl_carbon,) = carbonyl_neighbors
            n_alkyl_carbons = tuple(
                n.GetIdx() for n in neighbors if n.GetIdx() != carbonyl_carbon.GetIdx() and n.GetAtomicNum() == 6
            )
            amide_nitrogen_by_carbon.setdefault(carbonyl_carbon.GetIdx(), []).append(atom.GetIdx())
            n_alkyl_carbons_by_nitrogen[atom.GetIdx()] = n_alkyl_carbons
            n_hydroxy_by_nitrogen[atom.GetIdx()] = n_hydroxy_oxygen
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

    for carbon_idx in set(amide_oxygen_by_carbon) | set(amide_nitrogen_by_carbon):
        oxygens = amide_oxygen_by_carbon.get(carbon_idx, [])
        nitrogens = amide_nitrogen_by_carbon.get(carbon_idx, [])
        if len(oxygens) != 1 or len(nitrogens) != 1:
            raise UnsupportedStructure(
                "an oxygen/nitrogen pattern that isn't exactly one carbonyl "
                "oxygen and one primary-amide nitrogen on the same carbon "
                "is a more/less senior characteristic group than a plain "
                "primary amide (Table 3.3), which this module does not "
                "attempt to disambiguate"
            )
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "an amide carbon with more than one carbon neighbor is not "
                "a valid terminal amide carbon (a ketone-shaped carbon is "
                "out of scope for this module)"
            )
        amide_carbons.add(carbon_idx)

    if not amide_carbons:
        raise UnsupportedStructure(
            "no amide (-CON(R)(R')) group found; this module only "
            "handles amides"
        )
    if len(amide_carbons) > 1:
        raise UnsupportedStructure(
            "more than one amide group (a diamide) is out of scope for "
            "this module"
        )
    (amide_carbon,) = amide_carbons
    (amide_oxygen,) = amide_oxygen_by_carbon[amide_carbon]
    (amide_nitrogen,) = amide_nitrogen_by_carbon[amide_carbon]
    n_alkyl_carbons = n_alkyl_carbons_by_nitrogen[amide_nitrogen]
    n_hydroxy_oxygen = n_hydroxy_by_nitrogen[amide_nitrogen]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped, n_names=()):
    from ._amine import _add_n_names

    prefix = format_substituent_prefixes(_add_n_names(grouped, n_names))
    return prefix + name_from_substituents(chain_length, ene_locants, yne_locants, "amide")


def _candidate_key(chain_length, ene_locants, yne_locants, substituents, n_names=()):
    """Sort key implementing P-44.4.1.10 (ene/yne locants) ahead of P-45.2
    (substituent-prefix locants), most-preferred first. The amide group's
    own locant isn't part of this key: candidates are pre-filtered so the
    amide carbon always sits at C1 (see `_name_acyclic_amide`)."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped, n_names)
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




def _name_acyclic_amide(
    mol,
    amide_carbon,
    amide_nitrogen,
    excluded,
    n_alkyl_carbons,
    hydroxyls,
    bonds,
    stereo=None,
    extra_names=None,
    required_atoms=frozenset(),
    extra_excluded_carbons=frozenset(),
    phenyl_atoms=frozenset(),
    n_hydroxy_oxygen=None,
):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_carboxylic_acid.py`/`_aldehyde.py`'s identical treatment),
    and the winning candidate's own locants are used to format a
    "(<locant><R/S>,...)-" prefix onto the final name -- applied outermost,
    ahead of any N-alkyl prefix, since a stereodescriptor always sits at
    the very front of the complete name (P-91.3).

    `extra_names`/`required_atoms`: same coexisting-group injection point
    as `_ketone.py`'s `_name_acyclic_ketone`, reused by `_amide_amine.py`
    via `_coexisting_groups.py` -- both empty/None by default so existing
    callers are unaffected. `extra_excluded_carbons`: additional carbon
    atoms to remove from the principal-chain search graph before picking
    the longest chain, merged with the N-alkyl substituents' own excluded
    atoms below -- `_ether_amide.py` passes a coexisting ether's alkoxy-
    branch component here, the same reason and pattern as
    `_thiol.py`/`_ketone.py`/`_aldehyde.py`'s `carbon_graph` parameter (PR
    #428-430), just expressed as atoms-to-remove instead of a whole
    replacement graph since this function already builds its own."""
    graph = adjacency(mol)
    halogens = {
        **halogen_substituents(mol),
        **{o: "hydroxy" for o in hydroxyls},
        **(extra_names or {}),
    }
    full_carbon_graph = carbon_adjacency(mol)

    n_names = []
    n_substituent_atoms = set()
    for n_alkyl_c in n_alkyl_carbons:
        if n_alkyl_c in phenyl_atoms:
            n_names.append(name_branch(graph, n_alkyl_c, amide_nitrogen, {}, phenyl_atoms, mol=mol))
            n_substituent_atoms |= phenyl_atoms
            continue
        n_atoms, _ = bfs(full_carbon_graph, n_alkyl_c)
        n_atoms = set(n_atoms)
        if any(b[0] in n_atoms or b[1] in n_atoms for b in non_single_bonds(mol)):
            raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
        if any(next(iter(graph[o])) in n_atoms for o in hydroxyls):
            # A hydroxyl on the N-substituent is invisible to
            # `full_carbon_graph` (oxygen isn't a carbon), so it would
            # otherwise pass `name_branch` silently and get misnamed as a
            # plain, unsubstituted alkyl group -- the module docstring's
            # "plain, unsubstituted" N-substituent restriction is enforced
            # here explicitly (found via `specified_stereocenters`
            # correctly flagging this shape's stereocenters as partially
            # specified).
            raise UnsupportedStructure(
                "a substituted N-substituent (e.g. bearing a hydroxyl) is "
                "not supported yet; only a plain, unsubstituted alkyl "
                "N-substituent is in scope"
            )
        if any(nbr in n_atoms for h in halogen_substituents(mol) for nbr in graph[h]):
            # A halogen on the N-substituent is likewise invisible to
            # `full_carbon_graph` -- calling `name_branch` with an empty
            # halogens dict here would otherwise walk straight through it
            # as if it were a chain-extending atom (misnaming e.g.
            # '-CH2CH2Cl' as a 3-atom 'propyl' chain). Same
            # "plain, unsubstituted" restriction as the hydroxyl check
            # above, enforced explicitly rather than silently mishandled.
            raise UnsupportedStructure(
                "a substituted N-substituent (e.g. bearing a halogen) is "
                "not supported yet; only a plain, unsubstituted alkyl "
                "N-substituent is in scope"
            )
        n_names.append(name_branch(graph, n_alkyl_c, amide_nitrogen, {}, mol=mol))
        n_substituent_atoms |= n_atoms

    if n_hydroxy_oxygen is not None:
        # P-66.1.1.3.2: a hydroxamic acid is just an amide with a plain
        # 'hydroxy' N-substituent, cited via the same "N-" prefix
        # machinery below as any N-alkyl substituent.
        n_names.append(("hydroxy", False))

    # N-alkyl substituent carbons hang off the (excluded) amide nitrogen, not
    # off any acyl-chain carbon, so they form their own isolated component(s)
    # in the carbon-only graph; the whole subtree must be removed before
    # picking the longest chain, or a longer N-substituent (e.g.
    # N-butylacetamide) would be mistaken for the acyl chain itself.
    excluded_carbons = n_substituent_atoms | extra_excluded_carbons
    carbon_graph = {
        k: [n for n in v if n not in excluded_carbons]
        for k, v in full_carbon_graph.items()
        if k not in excluded_carbons
    }

    chains = all_chains(carbon_graph)
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if amide_carbon not in chain:
            continue
        if not required_atoms <= set(chain):
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            amide_carbon in c
            and required_atoms <= set(c)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the amide-bearing carbon (and/or a multiple bond) does not lie "
            "on a single longest carbon chain; a shorter principal chain is "
            "not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != amide_carbon:
                # The amide carbon must sit at C1 (see module docstring); a
                # direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents, n_names)
            if best_key is None or key < best_key:
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _ring_substituents(graph, ring_order, halogens, excluded, mol=None):
    """{ring position -> [substituent name, ...]}, mirroring
    `_nitrile.py`'s identically-named helper -- every branch hanging off
    a ring atom other than the amide carbon itself (in `excluded`) is a
    plain substituent prefix (alkyl/halogen)."""
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, amide_locant, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(ring_size, [], [], prefix, total_subs, "carboxamide", [amide_locant])


def _ring_candidate_key(ring_size, amide_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _ring_name_from_substituents(ring_size, amide_locant, grouped)
    return amide_locant, locant_set, citation_locants, name


def _name_ring_amide(mol, ring_atoms, stereo=None):
    """P-66.1.1.1.1.3: "The suffix 'carboxamide' is always used to name
    amides with the –CO-NH2 group attached to a ring" -- e.g.
    'cyclohexanecarboxamide'. Unlike the chain-parent 'amide' suffix
    above (where the -CONH2 carbon is always the parent's own C1), here
    the ring itself is the parent hydride and the -CONH2 carbon is a
    substituent atom hanging directly off one ring carbon, mirroring
    `_nitrile.py`'s `_name_ring_nitrile` construction.

    A saturated monocyclic all-carbon ring with exactly one N-unsubstituted
    -CONH2 hanging directly off one ring atom (no ring unsaturation, no
    standalone hydroxyl), plus any number of substituents (alkyl/halogen)
    on any ring atom, including the -CONH2-bearing atom itself.

    `stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- every stereocenter must lie on the ring
    itself (a stereocenter on the -CONH2 substituent branch is out of
    scope for now)."""
    amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen = (
        _validate_and_collect_amide(mol)
    )
    if n_alkyl_carbons:
        raise UnsupportedStructure("an N-alkyl-substituted ring amide is not supported yet")
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a ring amide is not "
            "supported yet"
        )

    graph = adjacency(mol)
    excluded_atoms = {amide_oxygen, amide_nitrogen}
    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_atoms and b[1] not in excluded_atoms and b[0] in ring_atoms and b[1] in ring_atoms
    ]
    if all_non_single:
        raise UnsupportedStructure(
            "an unsaturated ring alongside an amide substituent is not "
            "supported yet (see P-31.1.3)"
        )

    ring_neighbors = [n for n in graph[amide_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "an amide not directly attached to a single ring atom is not "
            "supported yet"
        )
    (ring_atom,) = ring_neighbors

    if stereo is not None and any(atom not in ring_atoms for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on the -CONH2 substituent branch rather than "
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
            amide_locant = position_of[ring_atom]
            substituents = _ring_substituents(graph, candidate, halogens, {amide_carbon}, mol=mol)
            key = _ring_candidate_key(ring_size, amide_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if n_hydroxy_oxygen is not None:
        separator = "-" if best_name[0].isdigit() else ""
        best_name = f"N-hydroxy{separator}{best_name}"

    if stereo:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _benzamide_name_from_substituents(grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if total_subs == 0:
        # 'benzamide' is a fully retained name (P-66.1.1.1.2.1) -- unlike
        # 'cyclohexanecarboxamide', there is no locant position to even
        # omit.
        return "benzamide"
    return f"{format_substituent_prefixes(grouped)}benzamide"


def _benzamide_candidate_key(amide_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzamide_name_from_substituents(grouped)
    return amide_locant, locant_set, citation_locants, name


def _name_benzamide(mol, ring_atoms, exempt_atoms=None):
    """P-66.1.1.1.2.1: 'benzamide' is one of only four retained amide
    names that are preferred IUPAC names and can be substituted -- itself
    the PIN for a -CONH2 hanging directly off one carbon of an otherwise-
    plain (or substituted) benzene ring -- e.g. 'benzamide' (PubChem CID
    239), '2-methylbenzamide' (CID 68140), '4-methylbenzamide' (CID
    69254). Structurally identical to `_name_ring_amide`'s saturated-ring
    case, mirroring `_nitrile.py`'s `_name_benzonitrile`, with the
    retained name 'benzamide' replacing 'cyclo' + alkane_name +
    'carboxamide' as the whole suffix unit (no locant is ever cited for
    the -CONH2 position itself)."""
    amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen = (
        _validate_and_collect_amide(mol, aromatic_ring_atoms=exempt_atoms or ring_atoms)
    )
    if n_alkyl_carbons or n_hydroxy_oxygen is not None:
        raise UnsupportedStructure("an N-substituted benzamide is not supported yet")
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside benzamide is not supported "
            "yet"
        )

    graph = adjacency(mol)
    excluded_atoms = {amide_oxygen, amide_nitrogen}
    ring_neighbors = [n for n in graph[amide_carbon] if n in ring_atoms]
    if len(ring_neighbors) != 1:
        raise UnsupportedStructure(
            "an amide not directly attached to a single ring atom is not "
            "supported yet"
        )
    (ring_atom,) = ring_neighbors
    other_ring_atom_branches = [n for n in graph[ring_atom] if n not in ring_atoms and n != amide_carbon]
    if other_ring_atom_branches:
        raise UnsupportedStructure(
            "a substituent on the same ring atom as the amide is not "
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
            amide_locant = position_of[ring_atom]
            substituents = _ring_substituents(graph, candidate, halogens, {amide_carbon}, mol=mol)
            key = _benzamide_candidate_key(amide_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    return best_name


def _name_phenyl_chain_amide(mol, ring_atoms):
    """Name a primary amide whose -CONH2 lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropanamide. The ring is
    cited as a 'phenyl' substituent prefix (via `name_branch`'s aromatic-
    ring recognition) on the chain, which is the parent hydride, mirroring
    `_aldehyde.py`'s `_name_phenyl_chain_aldehyde`. Narrower than the
    acyclic path above: no N-substitution, no coexisting standalone
    hydroxyl, no chain unsaturation, and no specified stereocenter -- each
    is a separate follow-up."""
    amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen = (
        _validate_and_collect_amide(mol, aromatic_ring_atoms=ring_atoms)
    )
    if n_alkyl_carbons or n_hydroxy_oxygen is not None:
        raise UnsupportedStructure(
            "an N-substituted amide alongside a benzene-ring substituent "
            "is not supported yet"
        )
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a benzene-ring-substituent "
            "amide chain is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "amide chain is not supported yet"
        )

    excluded = {amide_oxygen, amide_nitrogen}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded and b[1] not in excluded and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent amide "
            "chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_branch_attachments(mol, graph, rings)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain amide is not "
            "supported yet"
        )
    chain, branches = longest_branched_chain(graph, amide_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, [], [], grouped)


def _name_amide_with_n_phenyl(mol, ring_atoms):
    """Name a primary/secondary amide whose nitrogen carries a plain,
    unsubstituted benzene ring as a direct 'N-phenyl' substituent (e.g.
    acetanilide -> 'N-phenylacetamide', PubChem-confirmed), mirroring
    `_urea.py`'s/`_carbamate.py`'s identical N-phenyl pattern: the ring is
    just another N-substituent name (via `name_branch`'s aromatic-ring
    recognition), reusing `_name_acyclic_amide`'s existing N-prefix
    formatting unchanged. A second substituent sharing the same nitrogen
    (phenyl+alkyl, or two phenyls) is out of scope, same restriction as
    those modules; each is a separate follow-up."""
    amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen = (
        _validate_and_collect_amide(mol, aromatic_ring_atoms=ring_atoms)
    )
    if len(n_alkyl_carbons) != 1 or n_hydroxy_oxygen is not None:
        raise UnsupportedStructure(
            "a phenyl N-substituent alongside another substituent on the "
            "same nitrogen is not supported yet"
        )
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside an N-phenylamide is not "
            "supported yet"
        )
    stereo = specified_stereocenters(mol)
    excluded = {amide_oxygen, amide_nitrogen}
    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded and b[1] not in excluded and (b[0] not in ring_atoms or b[1] not in ring_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    return _name_acyclic_amide(
        mol,
        amide_carbon,
        amide_nitrogen,
        excluded,
        n_alkyl_carbons,
        hydroxyls,
        bonds,
        stereo,
        phenyl_atoms=ring_atoms,
    )


def name_amide(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        union = set().union(*aromatic_rings)
        amide_carbon, _, amide_nitrogen, *_ = _validate_and_collect_amide(mol, aromatic_ring_atoms=union)
        if any(n in union for n in adjacency(mol)[amide_nitrogen]):
            raise UnsupportedStructure("an N-aryl amide alongside a second aromatic ring is not supported yet")
        anchors = [amide_carbon]
        host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
        if host is not None:
            return _name_benzamide(mol, host, union)
        return _name_phenyl_chain_amide(mol, union)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            graph = adjacency(mol)
            attachment = ring_chain_attachment(graph, ring_atoms, set())
            if attachment is not None and mol.GetAtomWithIdx(attachment[1]).GetAtomicNum() == 7:
                return _name_amide_with_n_phenyl(mol, ring_atoms)
            amide_carbon, _, _, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen = _validate_and_collect_amide(
                mol, aromatic_ring_atoms=ring_atoms
            )
            if not n_alkyl_carbons and not hydroxyls and n_hydroxy_oxygen is None:
                graph = adjacency(mol)
                ring_neighbors = [n for n in graph[amide_carbon] if n in ring_atoms]
                if len(ring_neighbors) == 1:
                    return _name_benzamide(mol, ring_atoms)
            return _name_phenyl_chain_amide(mol, ring_atoms)
        amide_carbon, _, _, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen = _validate_and_collect_amide(mol)
        if not n_alkyl_carbons and not hydroxyls:
            graph = adjacency(mol)
            ring_neighbors = [n for n in graph[amide_carbon] if n in ring_atoms]
            if len(ring_neighbors) == 1:
                return _name_ring_amide(mol, ring_atoms, specified_stereocenters(mol))
    amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls, n_hydroxy_oxygen = (
        _validate_and_collect_amide(mol)
    )
    stereo = specified_stereocenters(mol)
    if carbon_on_ring(mol, adjacency(mol), [amide_carbon]):
        raise UnsupportedStructure(
            "this ring shape alongside an amide (a true lactam, N-alkyl "
            "substitution, more than one amide, a standalone hydroxyl, "
            "ring unsaturation, or a ring other than a single saturated "
            "monocyclic/benzene one) is out of scope for this module's "
            "'carboxamide' suffix path (P-66.1.1.1.1.3)"
        )
    excluded = {amide_oxygen, amide_nitrogen}
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
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

    return _name_acyclic_amide(
        mol,
        amide_carbon,
        amide_nitrogen,
        excluded,
        n_alkyl_carbons,
        hydroxyls,
        bonds,
        stereo,
        n_hydroxy_oxygen=n_hydroxy_oxygen,
    )
