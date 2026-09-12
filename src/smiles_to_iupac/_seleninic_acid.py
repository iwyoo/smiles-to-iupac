"""Naming of seleninic acids (the '-seleninic acid' suffix, -Se(=O)OH) on
acyclic saturated or unsaturated carbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-65.3.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf), the
  selenium analogue of the sulfinic acid entry in the same table (the same
  chalcogen-pair pattern as `_selenonic_acid.py`'s own relationship to
  `_sulfonic_acid.py`): 'seleninic acid' is the preselected suffix for
  -Se(=O)-OH, one oxidation state below selenonic acid -- the selenium
  carries one carbon, one double-bonded oxygen, and one hydroxyl oxygen
  (degree 3, not 4), structurally identical to `_sulfinic_acid.py`'s
  -S(=O)-OH with selenium instead of sulfur.
- Like 'sulfinic acid'/'selenonic acid', 'seleninic acid' is cited as the
  parent hydride name followed directly by the two-word suffix with no
  elision -- 'methane' + 'seleninic acid' -> 'methaneseleninic acid'
  (PubChem structure match).
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_sulfinic_acid.py`/`_selenonic_acid.py` apply, e.g. 'ethaneseleninic
  acid' (PubChem structure match).
- P-44.4.1.8 / P-45.2: the -Se(=O)OH locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as `_sulfinic_acid.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -Se(=O)OH suffix.
- P-92 stereocenters:
  like `_sulfinic_acid.py`'s sulfur, this module's seleninic selenium
  (-R, =O, -OH, a lone pair) is *itself* a potential stereocenter in
  essentially every real -Se(=O)OH molecule -- confirmed via RDKit's
  `Chem.FindPotentialStereo` on `CC(C)[Se](=O)O` (no chain stereocenter
  at all), which still flags the selenium atom. Same conclusion as
  sulfinic acid: this module only ever explicitly rejects a specified
  stereocenter (chain carbon or selenium alike) rather than attempting to
  cite one, for the same reasons documented in `_sulfinic_acid.py`'s
  module docstring.

Scope, deliberately narrow, mirroring `_selenonic_acid.py`'s own
chain-only first pass (no monocyclic seleninic acid has been found
registered on PubChem to verify that shape here): a single -Se(=O)OH on
an acyclic chain, with no other heteroatom anywhere in the molecule except
the seleninic acid group's own two oxygens. Explicitly out of scope (raise
`UnsupportedStructure`): any saturated ring, two or more -Se(=O)OH groups,
and a seleninic acid on a carbon that is also part of a C=C/C#C bond. One
*aromatic*-ring exception beyond the phenyl-chain path below:
`_name_benzeneseleninic_acid` names -Se(=O)OH directly on a benzene ring
carbon (with or without other ring substituents), e.g.
'benzeneseleninic acid', '4-methylbenzeneseleninic acid' (both PubChem
PUG REST matches), mirroring `_sulfonic_acid.py`'s
`_name_benzenesulfonic_acid` with the retained name 'benzene' as stem;
the -Se(=O)OH's own locant is never cited here.
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
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_SELENIUM = 34


def _seleninic_selenium_atoms(mol):
    """Selenium atoms shaped like a seleninic acid group: bonded to exactly
    one carbon, one double-bonded (terminal) oxygen, and one single-bonded
    hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _SELENIUM or atom.GetDegree() != 3:
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


def has_seleninic_acid_shape(mol) -> bool:
    return bool(_seleninic_selenium_atoms(mol))


def _validate_and_collect_seleninic_acids(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_seleninic_acid`'s benzene-
    ring-substituent path (see `_name_phenyl_chain_seleninic_acid`) can
    reuse this same validation for the rest of the molecule. Empty by
    default, so every other caller's behavior is unchanged."""
    selenium_atoms = _seleninic_selenium_atoms(mol)
    if not selenium_atoms:
        raise UnsupportedStructure(
            "no seleninic acid (-Se(=O)OH) group found; this module only "
            "handles seleninic acids"
        )
    if len(selenium_atoms) > 1:
        raise UnsupportedStructure(
            "more than one seleninic acid group is out of scope for this "
            "module"
        )
    seleninic_atom_idxs = set()
    for se in selenium_atoms:
        seleninic_atom_idxs.add(se.GetIdx())
        seleninic_atom_idxs.update(n.GetIdx() for n in se.GetNeighbors() if n.GetAtomicNum() == 8)

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
        elif atom.GetIdx() not in seleninic_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a seleninic acid group (P-65.3.1) "
                "and halogen substituents (P-35.2.1) are not supported "
                "yet -- in particular a coexisting carboxylic/sulfonic/"
                "sulfinic/selenonic acid or other characteristic group "
                "needs acid-vs-acid Table 3.3 seniority handling not yet "
                "implemented here"
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

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0 and not (
        ring_info.NumRings() == 1 and set(ring_info.AtomRings()[0]) == set(aromatic_ring_atoms)
    ):
        raise UnsupportedStructure(
            "rings are not supported yet, other than the separate benzene-"
            "ring-substituent path (this module otherwise only handles "
            "acyclic chains)"
        )

    (selenium,) = selenium_atoms
    (carbon,) = (n for n in selenium.GetNeighbors() if n.GetAtomicNum() == 6)
    return selenium.GetIdx(), carbon.GetIdx()


def _reject_eneseleninic_carbon(graph, seoh_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if seoh_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a seleninic acid on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
        )


def _suffix_body(ene_locants, yne_locants, seoh_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append(([seoh_locant], "seleninic acid"))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and elides_before(words[i + 1]):
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _name_from_substituents(chain_length, seoh_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "seleninic acid"

    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole -Se(=O)OH) in total omits the locant, e.g.
        # 'ethaneseleninic acid'.
        return alkane_name(2) + "seleninic acid"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, seoh_locant)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, seoh_locant, ene_locants, yne_locants, substituents):
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
    name = _name_from_substituents(chain_length, seoh_locant, ene_locants, yne_locants, grouped)
    return (
        (
            seoh_locant,
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


def _substituents_for_ring(graph, ring_order, halogens, excluded, mol=None):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _benzeneseleninic_acid_name_from_substituents(grouped):
    # Mirrors `_sulfonic_acid.py`'s `_benzenesulfonic_acid_name_from_substituents`:
    # the mancude ring's own numbering is always free to start at the
    # -Se(=O)OH carbon, so its locant is never cited even when other
    # substituents need theirs, e.g. '4-methylbenzeneseleninic acid'.
    if not grouped:
        return "benzeneseleninic acid"
    return f"{format_substituent_prefixes(grouped)}benzeneseleninic acid"


def _benzeneseleninic_acid_candidate_key(seoh_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _benzeneseleninic_acid_name_from_substituents(grouped)
    return seoh_locant, locant_set, citation_locants, name


def _name_benzeneseleninic_acid(mol, ring_atoms):
    """P-65.3.1: -Se(=O)OH attached directly to a benzene ring carbon --
    e.g. 'benzeneseleninic acid', '4-methylbenzeneseleninic acid' (both
    PubChem PUG REST matches). Mirrors `_sulfonic_acid.py`'s
    `_name_benzenesulfonic_acid` with the retained name 'benzene' as
    stem; the -Se(=O)OH's own locant is never cited here. A specified
    stereocenter is always rejected (module docstring: the seleninic
    selenium is itself a potential stereocenter with no established way
    to cite it)."""
    selenium_idx, seoh_carbon = _validate_and_collect_seleninic_acids(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter (chain carbon or the seleninic "
            "selenium itself) is not supported yet for seleninic acids "
            "(see P-92, module docstring)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {selenium_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            seoh_locant = position_of[seoh_carbon]
            substituents = _substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _benzeneseleninic_acid_candidate_key(seoh_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_seleninic_acid(mol, ring_atoms):
    """Name a seleninic acid whose -Se(=O)OH lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. phenylmethaneseleninic acid
    (PubChem CID 125616). The ring is cited as a 'phenyl' substituent
    prefix (via `name_branch`'s aromatic-ring recognition) on the chain,
    which is the parent hydride, mirroring `_sulfinic_acid.py`'s
    `_name_phenyl_chain_sulfinic_acid`. Narrower than the acyclic path
    above: no chain unsaturation."""
    selenium_idx, seoh_carbon = _validate_and_collect_seleninic_acids(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol) is not None:
        # See module docstring: the seleninic selenium is itself a
        # potential stereocenter in virtually every real -Se(=O)OH
        # molecule, and this project has no established way to cite one --
        # reject unconditionally, same as the acyclic path below.
        raise UnsupportedStructure(
            "a specified stereocenter (chain carbon or the seleninic "
            "selenium itself) is not supported yet for seleninic acids "
            "(see P-92, module docstring)"
        )
    excluded = {selenium_idx}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if selenium_idx not in (b[0], b[1]) and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "seleninic acid chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain seleninic acid is not supported yet"
        )
    chain, branches = longest_branched_chain_through(graph, seoh_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    halogens = halogen_substituents(mol)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        seoh_locant = position_of[seoh_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, seoh_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_seleninic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            _, seoh_carbon = _validate_and_collect_seleninic_acids(mol, aromatic_ring_atoms=ring_atoms)
            if seoh_carbon in ring_atoms:
                return _name_benzeneseleninic_acid(mol, ring_atoms)
            return _name_phenyl_chain_seleninic_acid(mol, ring_atoms)
    selenium_idx, seoh_carbon = _validate_and_collect_seleninic_acids(mol)
    if specified_stereocenters(mol) is not None:
        # Unlike `_sulfonic_acid.py`'s sulfur, this module's seleninic
        # selenium is itself a potential stereocenter in virtually every
        # real -Se(=O)OH molecule (module docstring), and this project has
        # no established way to cite a heteroatom-centered
        # stereodescriptor -- explicitly reject rather than silently drop
        # the marker (P-92), same as `_sulfinic_acid.py`.
        raise UnsupportedStructure(
            "a specified stereocenter (chain carbon or the seleninic "
            "selenium itself) is not supported yet for seleninic acids "
            "(see P-92, module docstring)"
        )
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and selenium_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 1:
        # The one Se=O double bond is always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_eneseleninic_carbon(graph, seoh_carbon, bonds)

    halogens = halogen_substituents(mol)
    excluded = {selenium_idx}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if seoh_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the seleninic-acid-bearing carbon (and/or a multiple bond) "
            "does not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            seoh_locant = position_of[seoh_carbon]
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, seoh_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
