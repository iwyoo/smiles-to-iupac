"""Naming of sulfinic acids (the '-sulfinic acid' suffix, -SO2H) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.3.1, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  'sulfinic acid' is the preselected suffix for -S(=O)-OH, one oxidation
  state below sulfonic acid (`_sulfonic_acid.py`) -- the sulfur carries one
  carbon, one double-bonded oxygen, and one hydroxyl oxygen (degree 3, not
  4). Table 3.3's acid cluster ranks sulfinic acid junior to sulfonic acid;
  this module doesn't implement that acid-vs-acid seniority (out of scope,
  see below -- this module only ever handles a single sulfinic acid group
  in isolation).
- Like 'sulfonic acid' (and unlike 'ol'/'thiol'/'amine'), 'sulfinic acid'
  is cited as the parent hydride name followed directly by the two-word
  suffix with no elision -- 'ethane' + 'sulfinic acid' -> 'ethanesulfinic
  acid'.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_sulfonic_acid.py` apply, e.g. 'ethanesulfinic acid'.
- P-44.4.1.8 / P-45.2: the -SO2H locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as every other suffix module here.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -SO2H suffix.
- P-92/P-93.3.4.1 stereocenters: unlike `_sulfonic_acid.py`'s sulfur (two
  identical =O, never stereogenic), this module's sulfinic sulfur
  (-R, =O, -OH, a lone pair) is *itself* a potential stereocenter in
  essentially every real -SO2H molecule -- confirmed via RDKit's
  `Chem.FindPotentialStereo` on e.g. plain `CC(C)S(=O)O` (no chain
  stereocenter at all), which still flags the sulfur atom. That means any
  input with a specified *chain* stereocenter almost always has this
  second, unspecified sulfur stereocenter riding along, and
  `_common.specified_stereocenters` correctly rejects that combination as
  partially specified (P-92) rather than silently ignoring either one.
  When the sulfur is the molecule's *sole* specified stereocenter,
  `_common.heteroatom_stereo_prefix` cites it with a bare "(R)-"/"(S)-"
  prefix -- P-93.3.4.1 assigns an ordinary R/S descriptor to this kind of
  trigonal pyramidal center, in the manner described for tetrahedral
  centers, e.g. "(R)-propane-1-sulfinic acid". A specified sulfur
  stereocenter *combined with* a specified chain/ring carbon stereocenter
  remains out of scope (two descriptors to combine into one citation
  group is a separate design question).

Scope, deliberately narrow (mirrors `_sulfonic_acid.py`'s own first pass):
a single -SO2H on an acyclic chain or on a single saturated carbon ring
(monocyclic, mirroring `_sulfonic_acid.py`'s own monocyclic support),
with no other heteroatom anywhere in the molecule except the sulfinic
acid group's own two oxygens -- acid-vs-acid seniority and any other
Table 3.3 seniority coexistence is future work, tracked under
`multi-carbonyl-seniority.md`.

P-31.1.3: a monocyclic ring bearing a sulfinic acid and exactly one C=C
ring double bond -- e.g. 'cyclohex-2-ene-1-sulfinic acid',
'cyclohex-3-ene-1-sulfinic acid', both confirmed via PubChem PUG REST.
Mirrors `_sulfonic_acid.py`'s identical extension; deliberately narrow: a
ring triple bond, and any other substituent alongside the ring double
bond, are both still explicitly rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): a specified
stereocenter on a von Baeyer polycyclic or spiro skeleton (a single
sulfinic acid on such a skeleton is supported, P-23/P-24 numbering
integration via `_polycyclic_suffix.py` with `elide_e=False`, mirroring
`_thiol.py`); unsaturation reaching outside the ring or a ring triple
bond, an -SO2H on a substituent branch off an otherwise-unsubstituted
*saturated* ring, two or more -SO2H groups, and a sulfinic acid on a
carbon that is also part of a C=C/C#C bond. Two *aromatic*-ring cases:
`_name_phenyl_chain_sulfinic_acid` names a -SO2H chain hanging off a
single plain, unsubstituted benzene ring (e.g.
'3-phenylpropane-1-sulfinic acid'), mirroring `_sulfonic_acid.py`'s
identical benzene-ring-substituent path -- narrower than the acyclic
path: no chain unsaturation. `_name_benzenesulfinic_acid` names -SO2H
directly on a benzene ring carbon (with or without other ring
substituents), e.g. 'benzenesulfinic acid' (PubChem CID 12057),
'2-methylbenzenesulfinic acid' (CID 12661295), mirroring
`_sulfonic_acid.py`'s `_name_benzenesulfonic_acid` with the retained
name 'benzene' as stem -- the -SO2H's own locant is never cited here,
unlike the cycloalkane case.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locant,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    heteroatom_stereo_prefix,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    name_from_substituents,
    non_single_bonds,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    ring_cycle,
    ring_name_from_substituents,
    substituent_locant_set_and_citation,
)
from ._bicyclic import find_bicyclic_core
from ._numerals import alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._substituents import (
    substituents_for_ring,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
    substituents_for_chain,
)

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _sulfinic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfinic acid group: bonded to exactly one
    carbon, one double-bonded (terminal) oxygen, and one single-bonded
    hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 2:
            continue
        double_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0]
        hydroxyl_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(double_os) != 1 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_sulfinic_acid_shape(mol) -> bool:
    return bool(_sulfinic_sulfur_atoms(mol))


def _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_sulfinic_acid`'s benzene-
    ring-substituent path (see `_name_phenyl_chain_sulfinic_acid`) can
    reuse this same validation for the rest of the molecule. Empty by
    default, so every other caller's behavior is unchanged."""
    sulfur_atoms = _sulfinic_sulfur_atoms(mol)
    if not sulfur_atoms:
        raise UnsupportedStructure(
            "no sulfinic acid (-SO2H) group found; this module only "
            "handles sulfinic acids"
        )
    if len(sulfur_atoms) > 1:
        raise UnsupportedStructure(
            "more than one sulfinic acid group is out of scope for this "
            "module"
        )
    sulfinic_atom_idxs = set()
    for s in sulfur_atoms:
        sulfinic_atom_idxs.add(s.GetIdx())
        sulfinic_atom_idxs.update(n.GetIdx() for n in s.GetNeighbors() if n.GetAtomicNum() == 8)

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
        elif atom.GetIdx() not in sulfinic_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a sulfinic acid group (P-65.3.1) "
                "and halogen substituents (P-35.2.1) are not supported "
                "yet -- in particular a coexisting carboxylic/sulfonic "
                "acid or other characteristic group needs acid-vs-acid "
                "Table 3.3 seniority handling not yet implemented here"
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
    return sulfur.GetIdx(), carbon.GetIdx()


def _reject_enesulfinic_carbon(graph, so2h_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if so2h_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a sulfinic acid on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
        )


def _name_from_substituents(chain_length, so2h_locant, ene_locants, yne_locants, grouped):
    # Only chain_length == 1 omits a substituent prefix's own locant too
    # (see `_alcohol.py`'s equivalent comment).
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    return prefix + name_from_substituents(chain_length, ene_locants, yne_locants, "sulfinic acid", [so2h_locant])


def _candidate_key(chain_length, so2h_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, so2h_locant, ene_locants, yne_locants, grouped)
    return (
        (
            so2h_locant,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )

def _ring_name_from_substituents(ring_size, so2h_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size, ene_locants, yne_locants, prefix, total_subs, "sulfinic acid", [so2h_locant]
    )


def _ring_candidate_key(ring_size, so2h_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, so2h_locant, ene_locants, yne_locants, grouped)
    return so2h_locant, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _name_cyclic_sulfinic_acid(mol, sulfur_idx, so2h_carbon, bonds=()):
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
            "sulfinic acid is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2h_locant = position_of[so2h_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, so2h_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _benzenesulfinic_acid_name_from_substituents(grouped):
    # Unlike the cycloalkane case, the mancude ring's own numbering is
    # always free to start at the -SO2H carbon (P-14.3.3-style), so its
    # locant is never cited even when other substituents need theirs,
    # e.g. '2-methylbenzenesulfinic acid' (PubChem CID 12661295), not
    # '2-methylbenzene-1-sulfinic acid' -- mirrors `_sulfonic_acid.py`'s
    # identical 'benzenesulfonic acid' treatment.
    if not grouped:
        return "benzenesulfinic acid"
    return f"{format_substituent_prefixes(grouped)}benzenesulfinic acid"


def _benzenesulfinic_acid_candidate_key(so2h_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzenesulfinic_acid_name_from_substituents(grouped)
    return so2h_locant, locant_set, citation_locants, name


def _name_benzenesulfinic_acid(mol, ring_atoms):
    """P-65.3.1: -SO2H attached directly to a benzene ring carbon -- e.g.
    'benzenesulfinic acid' (PubChem CID 12057), '2-methylbenzenesulfinic
    acid' (CID 12661295). Mirrors `_sulfonic_acid.py`'s
    `_name_benzenesulfonic_acid` exactly, with the retained name
    'benzene' as stem in place of 'cyclo' + alkane_name; an aromatic ring
    has no ene/yne locants of its own."""
    sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=ring_atoms)
    stereo_prefix = heteroatom_stereo_prefix(mol, sulfur_idx) or ""

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2h_locant = position_of[so2h_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _benzenesulfinic_acid_candidate_key(so2h_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return stereo_prefix + best_name


def _name_phenyl_chain_sulfinic_acid(mol, ring_atoms):
    """Name a sulfinic acid whose -SO2H lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropane-1-sulfinic acid.
    The ring is cited as a 'phenyl' substituent prefix (via
    `name_branch`'s aromatic-ring recognition) on the chain, which is the
    parent hydride, mirroring `_sulfonic_acid.py`'s
    `_name_phenyl_chain_sulfonic_acid`. Narrower than the acyclic path
    above: no chain unsaturation -- a separate follow-up (see
    tasks/phenyl-substituent-on-sulfinic-acid-chain.md's scope note)."""
    sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=ring_atoms)
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
            "sulfinic acid chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain sulfinic acid is "
            "not supported yet"
        )
    chain, branches = longest_branched_chain_through(graph, so2h_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        so2h_locant = position_of[so2h_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, so2h_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return stereo_prefix + best_name


def _name_ring_substituent_chain_sulfinic_acid(mol, sulfur_idx, so2h_carbon):
    """Name a sulfinic acid whose -SO2H lies entirely on a single branched
    chain hanging off one atom of an otherwise-plain saturated monocyclic
    ring (the ring itself bears no sulfinic acid) -- e.g.
    cyclohexylmethanesulfinic acid. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_sulfinic_acid` above and
    `_sulfonic_acid.py`'s `_name_ring_substituent_chain_sulfonic_acid`."""
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

    chain, branches = longest_branched_chain_through(graph, so2h_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
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
        so2h_locant = position_of[so2h_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, so2h_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_von_baeyer_or_spiro_sulfinic_acid(mol, sulfur_idx, so2h_carbon, bonds, stereo_prefix):
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single sulfinic acid (-S(=O)OH) suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`'s
    `elide_e=False` ('sulfinic acid' begins with a consonant, P-16.3.3,
    same as `_thiol.py`'s 'thiol' -- e.g. 'bicyclo[2.2.1]heptane-2-
    sulfinic acid', PubChem CID 20739363). Mirrors `_sulfinamide.py`'s
    own bicyclic/polycyclic-before-spiro dispatch order and restrictions:
    exactly one sulfinic acid on the ring system itself, no ring
    unsaturation, no specified stereocenter (on the ring or the sulfinic
    sulfur itself)."""
    if bonds:
        raise UnsupportedStructure(
            "an unsaturated von Baeyer bicyclic/polycyclic or monospiro "
            "ring system is not supported yet (see P-31.1.4/P-31.1.5)"
        )
    if stereo_prefix:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a von Baeyer bicyclic/"
            "polycyclic or monospiro sulfinic acid is not supported yet "
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
            so2h_carbon,
            {sulfur_idx},
            "sulfinic acid",
            "sulfinic acid",
            bicyclic_core,
            polycyclic_core,
            von_baeyer_ring_count,
            elide_e=False,
        )

    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return name_monospiro_suffix(
            mol, so2h_carbon, {sulfur_idx}, "sulfinic acid", "sulfinic acid", spiro_atom, elide_e=False
        )

    raise UnsupportedStructure(
        "polycyclic and fused-ring sulfinic acids are not supported yet "
        "(P-23/P-25 numbering integration with a suffix group is future "
        "work)"
    )


def name_sulfinic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol, aromatic_ring_atoms=ring_atoms)
            if so2h_carbon in ring_atoms:
                return _name_benzenesulfinic_acid(mol, ring_atoms)
            return _name_phenyl_chain_sulfinic_acid(mol, ring_atoms)
    sulfur_idx, so2h_carbon = _validate_and_collect_sulfinic_acids(mol)
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
    _reject_enesulfinic_carbon(graph, so2h_carbon, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        return _name_von_baeyer_or_spiro_sulfinic_acid(mol, sulfur_idx, so2h_carbon, bonds, stereo_prefix)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic sulfinic "
                "acid is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == _YNE_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a sulfinic "
                "acid is not supported yet -- only a ring double bond is "
                "in scope for this first pass (see P-31.1.3)"
            )
        if so2h_carbon not in ring_atoms:
            if not bonds:
                return stereo_prefix + _name_ring_substituent_chain_sulfinic_acid(mol, sulfur_idx, so2h_carbon)
            raise UnsupportedStructure(
                "a sulfinic acid on a substituent branch chain rather "
                "than the ring itself is not supported yet"
            )
        return stereo_prefix + _name_cyclic_sulfinic_acid(mol, sulfur_idx, so2h_carbon, bonds)

    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if so2h_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the sulfinic-acid-bearing carbon (and/or a multiple bond) "
            "does not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so2h_locant = position_of[so2h_carbon]
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, so2h_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return stereo_prefix + best_name
