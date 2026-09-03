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
  nothing, for methanamide/formamide, HCONH2). So the parent chain's
  numbering is never a locant-minimization choice for the amide group
  itself: whichever chain end carries the amide carbon simply becomes C1.
- P-14.3.3: the suffix locant for -CONH2 is never cited, since it is always
  fully determined by the amide carbon being a chain terminus with no other
  possible position, e.g. 'ethanamide' (not 'ethan-1-amide').
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
- The amide nitrogen may carry zero, one, or two plain, unbranched,
  unsubstituted, saturated alkyl substituents, each cited as its own
  'N-'-prefixed substituent directly ahead of the acyl stem, in alphabetical
  order, with a 'di' multiplying prefix (and a single shared 'N,N-' pair)
  when both are identical -- exactly `_carbamate.py`'s own N,N-disubstitution
  extension, reusing `alkyl_name` directly (not `name_branch`, see the P-29
  blocker) since both substituents are restricted to an unbranched chain.
  Confirmed via PubChem: 'N-methylacetamide' (CC(=O)NC), 'N,N-
  dimethylacetamide' (CC(=O)N(C)C), 'N-ethyl-N-methylacetamide'
  (CC(=O)N(C)CC).

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
- An N-substituent that is branched, unsaturated, or ring-bearing.
- An amide on/in a ring (a lactam) - a separate module's territory.
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
    bond_locant,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    linear_branch,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ordered_chain,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name, alkyl_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _is_carbonyl_carbon(mol, carbon_atom):
    return any(
        o.GetAtomicNum() == 8
        and o.GetDegree() == 1
        and mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        for o in carbon_atom.GetNeighbors()
    )


def has_amide_shape(mol) -> bool:
    """True if some carbon carries a doubly-bonded, monovalent carbonyl
    oxygen and a singly-bonded nitrogen with 0-2 carbon substituents, no
    other heavy-atom neighbor, and no *other* carbonyl-carbon neighbor (a
    -CON(R)(R') pattern, R/R' either H or an unbranched alkyl carbon),
    regardless of whether the rest of the molecule is in scope. Used by
    `core.py` to route ahead of the aldehyde/ketone dispatch, since an amide
    carbon would otherwise look aldehyde-shaped to those modules (both have
    exactly one carbon neighbor besides the carbonyl). A nitrogen bonded to
    two carbonyl carbons (a symmetric imide) is excluded here so `core.py`'s
    later `has_imide_shape` check still gets a chance at it."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if not _is_carbonyl_carbon(mol, atom):
            continue
        has_amide_n = any(
            n.GetAtomicNum() == 7
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and all(nn.GetAtomicNum() == 6 for nn in n.GetNeighbors() if nn.GetIdx() != atom.GetIdx())
            and sum(1 for nn in n.GetNeighbors() if _is_carbonyl_carbon(mol, nn)) == 1
            for n in atom.GetNeighbors()
        )
        if has_amide_n:
            return True
    return False


def _validate_and_collect_amide(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons,
    hydroxyls): the single -CON(R)(R') carbon/oxygen/nitrogen atom indices,
    a tuple of 0-2 N-alkyl substituent carbon indices, and the set of any
    coexisting standalone hydroxyl-oxygen atom indices.

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
            if any(n.GetAtomicNum() != 6 for n in neighbors):
                raise UnsupportedStructure(
                    "an amide nitrogen bonded to anything other than "
                    "carbon is out of scope for this module"
                )
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

            carbonyl_neighbors = [n for n in neighbors if _is_carbonyl_carbon(mol, n)]
            if len(carbonyl_neighbors) != 1:
                raise UnsupportedStructure(
                    "a nitrogen bonded to zero or multiple carbonyl carbons "
                    "is not a valid amide nitrogen for this module"
                )
            (carbonyl_carbon,) = carbonyl_neighbors
            n_alkyl_carbons = tuple(n.GetIdx() for n in neighbors if n.GetIdx() != carbonyl_carbon.GetIdx())
            amide_nitrogen_by_carbon.setdefault(carbonyl_carbon.GetIdx(), []).append(atom.GetIdx())
            n_alkyl_carbons_by_nitrogen[atom.GetIdx()] = n_alkyl_carbons
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
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls


def _suffix_body(ene_locants, yne_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'amide' ending
    (e.g. '2-enamide'); the amide group's own locant is never cited
    (P-14.3.3, see module docstring)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    words = [word for _, word in segments] + ["amide"]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    if segments:
        locant_parts = [
            f"{','.join(str(loc) for loc in locants)}-{word}"
            for (locants, _), word in zip(segments, words[:-1])
        ]
        body = "-".join(locant_parts) + words[-1]
    else:
        body = words[-1]
    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.10 (ene/yne locants) ahead of P-45.2
    (substituent-prefix locants), most-preferred first. The amide group's
    own locant isn't part of this key: candidates are pre-filtered so the
    amide carbon always sits at C1 (see `_name_acyclic_amide`)."""
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
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


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_amide(mol, amide_carbon, excluded, n_alkyl_carbons, hydroxyls, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_carboxylic_acid.py`/`_aldehyde.py`'s identical treatment),
    and the winning candidate's own locants are used to format a
    "(<locant><R/S>,...)-" prefix onto the final name -- applied outermost,
    ahead of any N-alkyl prefix, since a stereodescriptor always sits at
    the very front of the complete name (P-91.3)."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    full_carbon_graph = carbon_adjacency(mol)

    n_names = []
    n_substituent_atoms = set()
    for n_alkyl_c in n_alkyl_carbons:
        n_length = linear_branch(full_carbon_graph, n_alkyl_c, None)
        if n_length is None:
            raise UnsupportedStructure("a branched N-substituent is not supported yet")
        n_atoms = set()
        previous, current = None, n_alkyl_c
        while current is not None:
            n_atoms.add(current)
            neighbors = [n for n in full_carbon_graph[current] if n != previous]
            previous, current = current, (neighbors[0] if neighbors else None)
        if any(b[0] in n_atoms or b[1] in n_atoms for b in non_single_bonds(mol)):
            raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
        if any(next(iter(graph[o])) in n_atoms for o in hydroxyls):
            # A hydroxyl on the N-substituent is invisible to
            # `full_carbon_graph` (oxygen isn't a carbon), so it would
            # otherwise pass `linear_branch` silently and get misnamed as
            # a plain, unsubstituted alkyl group -- the module docstring's
            # "plain, unbranched, unsubstituted" N-substituent restriction
            # is enforced here explicitly (found via `specified_stereocenters`
            # correctly flagging this shape's stereocenters as partially
            # specified).
            raise UnsupportedStructure(
                "a substituted N-substituent (e.g. bearing a hydroxyl) is "
                "not supported yet; only a plain, unsubstituted alkyl "
                "N-substituent is in scope"
            )
        n_names.append(alkyl_name(n_length))
        n_substituent_atoms |= n_atoms

    # N-alkyl substituent carbons hang off the (excluded) amide nitrogen, not
    # off any acyl-chain carbon, so they form their own isolated component(s)
    # in the carbon-only graph; the whole subtree must be removed before
    # picking the longest chain, or a longer N-substituent (e.g.
    # N-butylacetamide) would be mistaken for the acyl chain itself.
    carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in n_substituent_atoms}

    chains = longest_chains(carbon_graph)
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if amide_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            amide_carbon in c and (not bonds or bond_locants(c, bonds) is not None) for c in chains
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
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != amide_carbon:
                # The amide carbon must sit at C1 (see module docstring); a
                # direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                best_key, best_name, best_position_of = key, name, position_of

    if n_alkyl_carbons:
        if len(n_names) == 2 and n_names[0] == n_names[1]:
            n_prefix = f"N,N-di{n_names[0]}"
        else:
            n_prefix = "-".join(f"N-{name}" for name in sorted(n_names))
        separator = "-" if best_name[0].isdigit() else ""
        best_name = f"{n_prefix}{separator}{best_name}"

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_amide(mol, ring_atoms):
    """Name a primary amide whose -CONH2 lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropanamide. The ring is
    cited as a 'phenyl' substituent prefix (via `name_branch`'s aromatic-
    ring recognition) on the chain, which is the parent hydride, mirroring
    `_aldehyde.py`'s `_name_phenyl_chain_aldehyde`. Narrower than the
    acyclic path above: no N-alkyl substitution, no coexisting standalone
    hydroxyl, no chain unsaturation, and no specified stereocenter -- each
    is a separate follow-up (see
    tasks/phenyl-substituent-on-amide-chain.md's scope note)."""
    amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls = _validate_and_collect_amide(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if n_alkyl_carbons:
        raise UnsupportedStructure(
            "an N-alkyl-substituted amide alongside a benzene-ring "
            "substituent is not supported yet"
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
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain amide is not supported yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, excluded)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside an "
            "amide is not supported yet"
        )
    if chain[-1] != amide_carbon:
        raise UnsupportedStructure(
            "the amide carbon must be the chain's far terminus from the "
            "benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an amide directly attached to the benzene ring (the "
            "'benzamide'-style naming) uses a separate construction, out "
            "of scope for this acyclic-chain-parent module"
        )

    ordered = list(reversed(chain))
    chain_length = len(ordered)
    position_of = {atom: i + 1 for i, atom in enumerate(ordered)}
    halogens = halogen_substituents(mol)
    substituents = {
        position_of[chain_root]: [name_branch(graph, ring_atom, chain_root, halogens, ring_atoms)]
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, [], [], grouped)


def name_amide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_amide(mol, ring_atoms)
    amide_carbon, amide_oxygen, amide_nitrogen, n_alkyl_carbons, hydroxyls = _validate_and_collect_amide(mol)
    stereo = specified_stereocenters(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an amide on/in a ring (a lactam) is out of scope for this "
            "acyclic-only module"
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

    return _name_acyclic_amide(mol, amide_carbon, excluded, n_alkyl_carbons, hydroxyls, bonds, stereo)
