"""Naming of diazonium cations (R-N#N+), per the IUPAC 2013
Recommendations ("the Blue Book"):

- Chapter P-7 (https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf, P-73): a
  diazonium cation is named by appending the suffix 'diazonium' directly
  to the parent hydride name, with no elision of the parent's terminal
  'e' (unlike the vowel-initial '-ol'/'-one' suffixes) -- 'methane' +
  'diazonium' -> 'methanediazonium'. This is structurally the same
  chain-search-and-locant machinery `_sulfonic_acid.py` already uses for
  its own '-SO3H' suffix, just with the diazonium nitrogen attached
  directly to the chain carbon instead of through an intervening sulfur.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_sulfonic_acid.py` apply (mononuclear parent, or a homogeneous
  two-carbon chain with exactly one substituent in total), e.g.
  'ethanediazonium'.
- Confirmed via PubChem structure match: `C[N+]#N` -> "methanediazonium",
  `CC[N+]#N` -> "ethanediazonium", `CCC[N+]#N` -> "propane-1-diazonium",
  `CC(C)[N+]#N` -> "propane-2-diazonium" (a branched chain, supported the
  same way `_sulfonic_acid.py` supports one), `ClCC[N+]#N` ->
  "2-chloroethanediazonium" (halogen coexistence).
- A single, otherwise-unsubstituted saturated monocyclic ring also works,
  mirroring `_sulfonic_acid.py`'s own monocyclic support (and
  `_radical.py`'s/`_carbenium.py`'s "cyclo" + parent name pattern):
  confirmed via PubChem structure match, `C1CCCCC1[N+]#N` ->
  "cyclohexanediazonium", `C1CCCC1[N+]#N` -> "cyclopentanediazonium". A
  ring-substituent case (a halogen or alkyl group elsewhere on the ring,
  e.g. `ClC1CCCCC1[N+]#N`) is NOT reachable here -- PubChem itself can't
  compute an IUPACName for that shape (an engine limitation, not a
  disproof), so it stays unverified and out of scope, unlike
  `_sulfonic_acid.py`'s own broader ring-substituent support.

-N#N+ directly on a benzene ring carbon (`_name_benzene_ring_diazonium`)
is also supported, with or without other ring substituents -- e.g.
'benzenediazonium' (PubChem CID 9718), '2-methylbenzenediazonium' (CID
192837). Mirrors `_carboxylic_acid.py`'s `_name_ring_carboxylic_acid`
ring-numbering search (the diazonium suffix locant is never cited, any
ring atom can be renumbered to position 1; other substituents get real
locants relative to that fixed reference), with the aromatic retained
name 'benzene' as stem instead of 'cyclo' + alkane_name. A chain-spacer
case (a -N#N+-bearing chain hanging off a plain benzene ring instead of
attaching directly) is a separate, already-supported construction --
see `_name_phenyl_chain_diazonium` below.

Scope, deliberately narrow, mirroring `_sulfonic_acid.py`'s own chain
scope: a single -N#N+ on an acyclic chain (branched, unbranched, or
unsaturated) with no other heteroatom anywhere in the molecule except the
diazonium group's own two nitrogens, OR a single, otherwise-unsubstituted
saturated monocyclic carbon ring, OR a single benzene ring with -N#N+
attached directly to one of its carbons (any number of other ring
substituents allowed). Explicitly out of scope (raise
`UnsupportedStructure`): a polycyclic/spiro/unsaturated non-benzene ring,
a substituent anywhere on an otherwise-unsubstituted saturated ring, two
or more diazonium groups, and a diazonium group on a carbon that is also
part of a C=C/C#C bond.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locant,
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
    ring_chain_attachment,
    ring_cycle,
    substituent_locant_set_and_citation,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _diazonium_nitrogens(mol):
    """Terminal, +1-charged nitrogens shaped like a diazonium group: triple
    bonded to a second, neutral, degree-1 nitrogen, and singly bonded to
    exactly one carbon."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1 or atom.GetDegree() != 2:
            continue
        if atom.GetIsotope() != 0:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(nitrogens) != 1:
            continue
        (carbon,) = carbons
        (terminal_n,) = nitrogens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), carbon.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), terminal_n.GetIdx()).GetBondTypeAsDouble() != _YNE_ORDER:
            continue
        if terminal_n.GetDegree() != 1 or terminal_n.GetFormalCharge() != 0 or terminal_n.GetIsotope() != 0:
            continue
        matches.append(atom)
    return matches


def has_diazonium_shape(mol) -> bool:
    return bool(_diazonium_nitrogens(mol))


def _validate_and_collect_diazonium(mol, aromatic_ring_atoms=frozenset()):
    """(diazonium_carbon_idx, diazonium_atom_idxs) after checking the
    molecule fits this module's scope (see module docstring).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection and the unconditional ring rejection below so
    `name_diazonium`'s benzene-ring-substituent path (see
    `_name_phenyl_chain_diazonium`) can reuse this same validation for the
    rest of the molecule. Empty by default, so every other caller's
    behavior is unchanged. Mirrors `_sulfonic_acid.py`'s equivalent
    aromatic-exemption pattern."""
    nitrogens = _diazonium_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure(
            "no diazonium (-N#N+) group found; this module only handles "
            "diazonium cations"
        )
    if len(nitrogens) > 1:
        raise UnsupportedStructure(
            "more than one diazonium group is out of scope for this module"
        )
    (nitrogen,) = nitrogens
    diazonium_atom_idxs = {nitrogen.GetIdx()}
    (terminal_n,) = (n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 7)
    diazonium_atom_idxs.add(terminal_n.GetIdx())

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
        elif atom.GetIdx() not in diazonium_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a diazonium group (P-73) and "
                "halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetIdx() not in diazonium_atom_idxs and (atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0):
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if not aromatic_ring_atoms and mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring other than a single, otherwise-unsubstituted saturated "
            "monocyclic carbon ring is not supported yet"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (carbon,) = (n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 6)
    return carbon.GetIdx(), diazonium_atom_idxs


def _unsubstituted_monocyclic_ring_diazonium_name(mol):
    """The name for a single -N#N+ on a single, otherwise-unsubstituted
    saturated monocyclic carbon ring (e.g. 'cyclohexanediazonium'), or
    None if the molecule isn't shaped like that at all."""
    nitrogens = _diazonium_nitrogens(mol)
    if len(nitrogens) != 1:
        return None
    (nitrogen,) = nitrogens
    (carbon,) = (n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 6)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    (ring_atoms,) = ring_info.AtomRings()
    if carbon.GetIdx() not in ring_atoms:
        return None
    if len(ring_atoms) != mol.GetNumAtoms() - 2:
        # every non-ring atom must be one of the diazonium group's own two
        # nitrogens -- no halogen or other substituent anywhere.
        return None

    for atom in mol.GetAtoms():
        if atom.GetIdx() in ring_atoms:
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None
            expected_degree = 3 if atom.GetIdx() == carbon.GetIdx() else 2
            if atom.GetDegree() != expected_degree:
                return None
        elif atom.GetIdx() not in (nitrogen.GetIdx(), next(n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 7)):
            return None
    diazonium_idxs = {nitrogen.GetIdx(), next(n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 7)}
    if any(b[0] not in diazonium_idxs or b[1] not in diazonium_idxs for b in non_single_bonds(mol)):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    return "cyclo" + alkane_name(len(ring_atoms)) + "diazonium"


def _reject_enediazonium_carbon(graph, diazonium_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if diazonium_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a diazonium group on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
        )


def _suffix_body(ene_locants, yne_locants, diazonium_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append(([diazonium_locant], "diazonium"))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and elides_before(words[i + 1]):
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _name_from_substituents(chain_length, diazonium_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "diazonium"

    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole diazonium group) in total omits the
        # locant, e.g. 'ethanediazonium'.
        return alkane_name(2) + "diazonium"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, diazonium_locant)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, diazonium_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, diazonium_locant, ene_locants, yne_locants, grouped)
    return (
        (
            diazonium_locant,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _substituents_for_chain(graph, chain, halogens, excluded, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _ring_substituents_diazonium(graph, ring_order, halogens, excluded, mol=None):
    """Mirrors `_sulfonic_acid.py`'s/`_carboxylic_acid.py`'s identically-
    named ring-substituent helpers: every branch hanging off a ring atom
    other than the diazonium group's own two nitrogens (in `excluded`) is
    a plain substituent prefix (alkyl/halogen)."""
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _candidate_key_diazonium(diazonium_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    if not grouped:
        name = "benzenediazonium"
    else:
        name = f"{format_substituent_prefixes(grouped)}benzenediazonium"
    return (diazonium_locant, locant_set, citation_locants, name), name


def _name_benzene_ring_diazonium(mol, ring_atoms):
    """P-73/P-14.3.3: -N#N+ attached directly to a benzene ring carbon --
    e.g. 'benzenediazonium' (PubChem CID 9718),
    '2-methylbenzenediazonium' (CID 192837). The diazonium suffix locant
    is never cited (any ring atom can be renumbered to position 1),
    mirroring `_carboxylic_acid.py`'s `_name_ring_carboxylic_acid`
    ring-numbering search, with the aromatic retained name 'benzene' as
    stem in place of 'cyclo' + alkane_name. Unlike that sp3-ring case, an
    aromatic ring carbon has no room for a second exocyclic substituent
    alongside the diazonium nitrogen (valence), so there is no
    same-ring-atom-substituent case to guard against here."""
    diazonium_carbon, excluded = _validate_and_collect_diazonium(mol, aromatic_ring_atoms=ring_atoms)
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
            diazonium_locant = position_of[diazonium_carbon]
            substituents = _ring_substituents_diazonium(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key_diazonium(diazonium_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def _name_phenyl_chain_diazonium(mol, ring_atoms):
    """Name a diazonium cation whose -N#N+ lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropane-1-diazonium. The
    ring is cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_sulfonic_acid.py`'s `_name_phenyl_chain_sulfonic_acid`.
    Narrower than the acyclic path above: no chain unsaturation."""
    diazonium_carbon, excluded = _validate_and_collect_diazonium(mol, aromatic_ring_atoms=ring_atoms)
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded and b[1] not in excluded and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "diazonium chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain diazonium is not supported yet"
        )
    chain, branches = longest_branched_chain_through(graph, diazonium_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    halogens = halogen_substituents(mol)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        diazonium_locant = position_of[diazonium_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, diazonium_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_diazonium(mol) -> str:
    ring_name = _unsubstituted_monocyclic_ring_diazonium_name(mol)
    if ring_name is not None:
        return ring_name

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            nitrogens = _diazonium_nitrogens(mol)
            if len(nitrogens) == 1:
                (nitrogen,) = nitrogens
                (diazonium_carbon,) = (
                    n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 6
                )
                if diazonium_carbon in ring_atoms:
                    return _name_benzene_ring_diazonium(mol, ring_atoms)
            return _name_phenyl_chain_diazonium(mol, ring_atoms)

    diazonium_carbon, excluded = _validate_and_collect_diazonium(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    # The diazonium N#N triple bond is the one non-single bond with both
    # atoms in `excluded`; anything else non-single must be a chain
    # ene/yne bond.
    bonds = [
        b for b in all_non_single
        if b[2] in (_ENE_ORDER, _YNE_ORDER) and b[0] not in excluded and b[1] not in excluded
    ]
    diazonium_bonds = [b for b in all_non_single if b[0] in excluded and b[1] in excluded]
    if len(bonds) + len(diazonium_bonds) != len(all_non_single) or len(diazonium_bonds) != 1:
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enediazonium_carbon(graph, diazonium_carbon, bonds)

    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if diazonium_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the diazonium-bearing carbon (and/or a multiple bond) does "
            "not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            diazonium_locant = position_of[diazonium_carbon]
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, diazonium_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
