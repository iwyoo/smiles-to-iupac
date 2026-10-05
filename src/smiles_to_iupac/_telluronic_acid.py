"""Naming of telluronic acids (the '-telluronic acid' suffix, -Te(=O)(=O)OH)
on acyclic saturated or unsaturated carbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-65.3.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf), the
  tellurium analogue of the sulfonic acid entry in the same table (the same
  chalcogen-quartet pattern already seen repeatedly in this project, e.g.
  the thioic/selenoic/telluroic acid and isothiocyanate/isoselenocyanate/
  isotellurocyanate modules): 'telluronic acid' is the preselected suffix
  for -Te(=O)(=O)-OH, structurally identical to `_sulfonic_acid.py`'s
  -S(=O)(=O)-OH (and `_selenonic_acid.py`'s -Se(=O)(=O)-OH) with tellurium
  instead.
- Like 'sulfonic acid'/'selenonic acid', 'telluronic acid' is cited as the
  parent hydride name followed directly by the two-word suffix with no
  elision -- 'methane' + 'telluronic acid' -> 'methanetelluronic acid'
  (PubChem structure match -- the only telluronic acid PubChem has
  registered at all; ethane/propane/branched-chain candidates all came
  back as CID 0). P-14.3.4.2(a)/(b)'s chain locant-omission rules
  themselves are NOT independently re-verified for tellurium here -- they
  are inherited, unchanged, from the identical mechanism already
  independently confirmed six times over in this project (`_thiol.py`,
  `_alcohol.py`, `_sulfonic_acid.py`, `_selenonic_acid.py`,
  `_sulfinic_acid.py`, `_seleninic_acid.py`), the same reduced-evidence
  bar `_selenol.py` already used for its own halogen-coexistence claim
  ("reviewed (eyeballed), not independently verified").
- P-44.4.1.8 / P-45.2: the -Te(=O)(=O)OH locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as `_sulfonic_acid.py`/`_selenonic_acid.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -Te(=O)(=O)OH suffix.

Scope, deliberately narrow, mirroring `_selenonic_acid.py`'s own
chain-only scope: a single -Te(=O)(=O)OH on an acyclic chain, with no
other heteroatom anywhere in the molecule except the telluronic acid
group's own three oxygens. Explicitly out of scope (raise
`UnsupportedStructure`): any saturated ring, two or more -Te(=O)(=O)OH
groups, and a telluronic acid on a carbon that is also part of a C=C/C#C
bond. One *aromatic*-ring exception: `_name_benzenetelluronic_acid` names
-Te(=O)(=O)OH directly on a benzene ring carbon (with or without other
ring substituents), e.g. 'benzenetelluronic acid',
'4-methylbenzenetelluronic acid' (both PubChem PUG REST IUPACName
matches), mirroring `_selenonic_acid.py`'s `_name_benzeneselenonic_acid`
with the retained name 'benzene' as stem; the -Te(=O)(=O)OH's own locant
is never cited here. A benzene-ring-substituent *chain* (mirroring
`_seleninic_acid.py`'s `_name_phenyl_chain_seleninic_acid`) is not
registered on PubChem for this suffix and stays out of scope, same as
`_selenonic_acid.py`.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    all_chains,
    bond_locant,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    lowest_locant_set,
    most_multiple_bonds,
    name_from_substituents,
    non_single_bonds,
    ring_cycle,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._substituents import format_substituent_prefixes, name_branch, ring_branch_stereo_display, substituents_for_ring, substituents_for_chain

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_TELLURIUM = 52


def _telluronic_tellurium_atoms(mol):
    """Tellurium atoms shaped like a telluronic acid group: bonded to exactly
    one carbon, two double-bonded (terminal) oxygens, and one single-bonded
    hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _TELLURIUM or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 3:
            continue
        double_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0]
        hydroxyl_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(double_os) != 2 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_telluronic_acid_shape(mol) -> bool:
    return bool(_telluronic_tellurium_atoms(mol))


def _validate_and_collect_telluronic_acids(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring -- exempted from the aromatic-atom rejection below so
    `name_telluronic_acid`'s benzene-ring-direct path (see
    `_name_benzenetelluronic_acid`) can reuse this same validation for the
    rest of the molecule. Empty by default, so every other caller's
    behavior is unchanged."""
    tellurium_atoms = _telluronic_tellurium_atoms(mol)
    if not tellurium_atoms:
        raise UnsupportedStructure(
            "no telluronic acid (-Te(=O)(=O)OH) group found; this module "
            "only handles telluronic acids"
        )
    if len(tellurium_atoms) > 1:
        raise UnsupportedStructure(
            "more than one telluronic acid group is out of scope for this "
            "module"
        )
    telluronic_atom_idxs = set()
    for te in tellurium_atoms:
        telluronic_atom_idxs.add(te.GetIdx())
        telluronic_atom_idxs.update(n.GetIdx() for n in te.GetNeighbors() if n.GetAtomicNum() == 8)

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
        elif atom.GetIdx() not in telluronic_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a telluronic acid group (P-65.3.1) "
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

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0 and not (
        ring_info.NumRings() == 1 and set(ring_info.AtomRings()[0]) == set(aromatic_ring_atoms)
    ):
        raise UnsupportedStructure(
            "rings other than a single plain benzene ring bearing the "
            "telluronic acid directly are not supported yet"
        )

    (tellurium,) = tellurium_atoms
    (carbon,) = (n for n in tellurium.GetNeighbors() if n.GetAtomicNum() == 6)
    return tellurium.GetIdx(), carbon.GetIdx()


def _reject_enetelluronic_carbon(graph, teo3h_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if teo3h_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a telluronic acid on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
        )


def _name_from_substituents(chain_length, teo3h_locant, ene_locants, yne_locants, grouped):
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    return prefix + name_from_substituents(chain_length, ene_locants, yne_locants, "telluronic acid", [teo3h_locant], substituted=bool(grouped))


def _candidate_key(chain_length, teo3h_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, teo3h_locant, ene_locants, yne_locants, grouped)
    return (
        (
            teo3h_locant,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )

def _benzenetelluronic_acid_name_from_substituents(grouped):
    # Mirrors `_selenonic_acid.py`'s
    # `_benzeneselenonic_acid_name_from_substituents`: the mancude ring's
    # own numbering is always free to start at the -Te(=O)(=O)OH carbon,
    # so its locant is never cited even when other substituents need
    # theirs, e.g. '4-methylbenzenetelluronic acid'.
    if not grouped:
        return "benzenetelluronic acid"
    return f"{format_substituent_prefixes(grouped)}benzenetelluronic acid"


def _benzenetelluronic_acid_candidate_key(teo3h_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzenetelluronic_acid_name_from_substituents(grouped)
    return teo3h_locant, locant_set, citation_locants, name


def _name_benzenetelluronic_acid(mol, ring_atoms):
    """P-65.3.1: -Te(=O)(=O)OH attached directly to a benzene ring carbon
    -- e.g. 'benzenetelluronic acid', '4-methylbenzenetelluronic acid'
    (both PubChem PUG REST IUPACName matches). Mirrors
    `_selenonic_acid.py`'s `_name_benzeneselenonic_acid` with the
    retained name 'benzene' as stem; the -Te(=O)(=O)OH's own locant is
    never cited here."""
    tellurium_idx, teo3h_carbon = _validate_and_collect_telluronic_acids(mol, aromatic_ring_atoms=ring_atoms)
    stereo = specified_stereocenters(mol)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {tellurium_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    branch_stereo = None
    if stereo:
        branch_stereo = ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a specified stereocenter alongside benzenetelluronic acid is "
                "not supported yet"
            )

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            teo3h_locant = position_of[teo3h_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            key = _benzenetelluronic_acid_candidate_key(teo3h_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_telluronic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzenetelluronic_acid(mol, ring_atoms)
    tellurium_idx, teo3h_carbon = _validate_and_collect_telluronic_acids(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and tellurium_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 2:
        # The two Se=O double bonds are always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enetelluronic_carbon(graph, teo3h_carbon, bonds)

    halogens = halogen_substituents(mol)
    excluded = {tellurium_idx}
    chains = all_chains(carbon_adjacency(mol))
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if teo3h_carbon not in chain:
            continue
        if stereo is not None and any(atom not in chain for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            teo3h_carbon in c for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the telluronic-acid-bearing carbon (and/or a multiple bond) "
            "does not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            teo3h_locant = position_of[teo3h_carbon]
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, teo3h_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
