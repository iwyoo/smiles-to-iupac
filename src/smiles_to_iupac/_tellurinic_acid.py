"""Naming of tellurinic acids (the '-tellurinic acid' suffix, -Te(=O)OH) on
acyclic saturated or unsaturated carbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-65.3.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf), the
  tellurium analogue of the sulfinic acid entry in the same table (the same
  chalcogen-pair pattern as `_seleninic_acid.py`'s own relationship to
  `_sulfinic_acid.py`): 'tellurinic acid' is the preselected suffix for
  -Te(=O)-OH, one oxidation state below telluronic acid (`_telluronic_acid.py`)
  -- the tellurium carries one carbon, one double-bonded oxygen, and one
  hydroxyl oxygen (degree 3, not 4), structurally identical to
  `_sulfinic_acid.py`'s -S(=O)-OH with tellurium instead of sulfur.
- Like 'sulfinic acid'/'seleninic acid', 'tellurinic acid' is cited as the
  parent hydride name followed directly by the two-word suffix with no
  elision -- 'methane' + 'tellurinic acid' -> 'methanetellurinic acid'
  (PubChem structure match -- the only tellurinic acid PubChem has
  registered at all; ethane/propane candidates come back as CID 0, the
  same sparse-data gap `_telluronic_acid.py` already documented). The
  chain locant-omission rules (P-14.3.4.2(a)/(b)) themselves are NOT
  independently re-verified for tellurium here -- inherited, unchanged,
  from the identical mechanism already independently confirmed by
  `_thiol.py`/`_alcohol.py`/`_sulfonic_acid.py`/`_selenonic_acid.py`/
  `_sulfinic_acid.py`/`_seleninic_acid.py`/`_telluronic_acid.py`, the same
  reduced-evidence bar `_selenol.py` first used for its own
  halogen-coexistence claim.
- P-44.4.1.8 / P-45.2: the -Te(=O)OH locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as `_sulfinic_acid.py`/`_telluronic_acid.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -Te(=O)OH suffix.
- P-91.3/P-92:
  unlike `_sulfinic_acid.py`'s sulfur or `_seleninic_acid.py`'s selenium,
  RDKit's `Chem.FindPotentialStereo` does not flag this module's
  tellurinic tellurium as a potential stereocenter at all, on any tried
  input -- this project delegates all CIP/stereocenter determination to
  RDKit, so per that same policy, a chain carbon stereocenter here is
  treated as the molecule's only stereo element, exactly as in
  `_sulfonic_acid.py`. A molecule with one or more *specified* tetrahedral
  stereocenters -- every one on the principal chain itself, no unspecified
  one alongside them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
  `_carboxylic_acid.py`/`_sulfonic_acid.py`.

Scope, deliberately narrow, mirroring `_telluronic_acid.py`'s own
chain-only scope: a single -Te(=O)OH on an acyclic chain, with no other
heteroatom anywhere in the molecule except the tellurinic acid group's
own two oxygens. Explicitly out of scope (raise `UnsupportedStructure`):
any saturated ring, two or more -Te(=O)OH groups, and a tellurinic acid
on a carbon that is also part of a C=C/C#C bond. One *aromatic*-ring
exception: `_name_benzenetellurinic_acid` names -Te(=O)OH directly on a
benzene ring carbon (with or without other ring substituents), e.g.
'benzenetellurinic acid', '4-methylbenzenetellurinic acid' (both PubChem
PUG REST IUPACName matches), mirroring `_telluronic_acid.py`'s
`_name_benzenetelluronic_acid`/`_selenonic_acid.py`'s
`_name_benzeneselenonic_acid` with the retained name 'benzene' as stem;
the -Te(=O)OH's own locant is never cited here. A benzene-ring-substituent
*chain* (mirroring `_seleninic_acid.py`'s
`_name_phenyl_chain_seleninic_acid`) is not registered on PubChem for
this suffix (CID 0 for phenylmethanetellurinic acid) and stays out of
scope.
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
    is_plain_benzene_ring,
    longest_chains,
    lowest_locant_set,
    non_single_bonds,
    ring_cycle,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    suffix_body,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_TELLURIUM = 52


def _tellurinic_tellurium_atoms(mol):
    """Tellurium atoms shaped like a tellurinic acid group: bonded to exactly
    one carbon, one double-bonded (terminal) oxygen, and one single-bonded
    hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _TELLURIUM or atom.GetDegree() != 3:
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


def has_tellurinic_acid_shape(mol) -> bool:
    return bool(_tellurinic_tellurium_atoms(mol))


def _validate_and_collect_tellurinic_acids(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring -- exempted from the aromatic-atom rejection below so
    `name_tellurinic_acid`'s benzene-ring-direct path (see
    `_name_benzenetellurinic_acid`) can reuse this same validation for the
    rest of the molecule. Empty by default, so every other caller's
    behavior is unchanged."""
    tellurium_atoms = _tellurinic_tellurium_atoms(mol)
    if not tellurium_atoms:
        raise UnsupportedStructure(
            "no tellurinic acid (-Te(=O)OH) group found; this module only "
            "handles tellurinic acids"
        )
    if len(tellurium_atoms) > 1:
        raise UnsupportedStructure(
            "more than one tellurinic acid group is out of scope for this "
            "module"
        )
    tellurinic_atom_idxs = set()
    for te in tellurium_atoms:
        tellurinic_atom_idxs.add(te.GetIdx())
        tellurinic_atom_idxs.update(n.GetIdx() for n in te.GetNeighbors() if n.GetAtomicNum() == 8)

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
        elif atom.GetIdx() not in tellurinic_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a tellurinic acid group (P-65.3.1) "
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
            "rings other than a single plain benzene ring bearing the "
            "tellurinic acid directly are not supported yet"
        )

    (tellurium,) = tellurium_atoms
    (carbon,) = (n for n in tellurium.GetNeighbors() if n.GetAtomicNum() == 6)
    return tellurium.GetIdx(), carbon.GetIdx()


def _reject_enetellurinic_carbon(graph, teoh_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if teoh_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a tellurinic acid on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
        )


def _name_from_substituents(chain_length, teoh_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "tellurinic acid"

    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole -Te(=O)OH) in total omits the locant, e.g.
        # 'ethanetellurinic acid'.
        return alkane_name(2) + "tellurinic acid"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = suffix_body(ene_locants, yne_locants, "tellurinic acid", [teoh_locant])[0]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, teoh_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, teoh_locant, ene_locants, yne_locants, grouped)
    return (
        (
            teoh_locant,
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


def _benzenetellurinic_acid_name_from_substituents(grouped):
    # Mirrors `_selenonic_acid.py`'s
    # `_benzeneselenonic_acid_name_from_substituents`: the mancude ring's
    # own numbering is always free to start at the -Te(=O)OH carbon, so
    # its locant is never cited even when other substituents need theirs,
    # e.g. '4-methylbenzenetellurinic acid'.
    if not grouped:
        return "benzenetellurinic acid"
    return f"{format_substituent_prefixes(grouped)}benzenetellurinic acid"


def _benzenetellurinic_acid_candidate_key(teoh_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzenetellurinic_acid_name_from_substituents(grouped)
    return teoh_locant, locant_set, citation_locants, name


def _name_benzenetellurinic_acid(mol, ring_atoms):
    """P-65.3.1: -Te(=O)OH attached directly to a benzene ring carbon --
    e.g. 'benzenetellurinic acid', '4-methylbenzenetellurinic acid' (both
    PubChem PUG REST IUPACName matches). Mirrors `_telluronic_acid.py`'s
    `_name_benzenetelluronic_acid` with the retained name 'benzene' as
    stem; the -Te(=O)OH's own locant is never cited here."""
    tellurium_idx, teoh_carbon = _validate_and_collect_tellurinic_acids(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter alongside benzenetellurinic acid is "
            "not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {tellurium_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            teoh_locant = position_of[teoh_carbon]
            substituents = _substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _benzenetellurinic_acid_candidate_key(teoh_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_tellurinic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzenetellurinic_acid(mol, ring_atoms)
    tellurium_idx, teoh_carbon = _validate_and_collect_tellurinic_acids(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and tellurium_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 1:
        # The one Se=O double bond is always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enetellurinic_carbon(graph, teoh_carbon, bonds)

    halogens = halogen_substituents(mol)
    excluded = {tellurium_idx}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if teoh_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        if stereo is not None and any(atom not in chain for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            teoh_carbon in c and (not bonds or bond_locants(c, bonds) is not None) for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the tellurinic-acid-bearing carbon (and/or a multiple bond) "
            "does not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            teoh_locant = position_of[teoh_carbon]
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, teoh_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
