"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
one or more unsubstituted sulfonamide (-SO2NH2) groups on the same
acyclic saturated carbon chain, per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-41/P-43 (`_seniority.py`): carboxylic acid (Table 4.4 class 1)
  outranks sulfonamide (class 19), so a coexisting unsubstituted
  sulfonamide is demoted to the 'sulfamoyl' substituent prefix
  (P-65.3.1) instead of its own '-sulfonamide' suffix -- e.g.
  'OC(=O)CS(=O)(=O)N' -> '2-sulfamoylacetic acid' (PubChem CID 11083961
  confirms the systematic '2-sulfamoylacetic acid' form; this module
  always uses the systematic '...anoic acid' stem rather than a retained
  name, mirroring `_carboxylic_acid.py`'s and
  `_carboxylic_acid_sulfinic_acid.py`'s own established convention once
  any substituent is present). Combines the parent-chain pattern from
  `_carboxylic_acid_sulfinic_acid.py` (PR #315, carboxylic acid wins the
  suffix) with the *unsubstituted-sulfonamide-only* shape check from
  `_sulfonic_acid_sulfonamide.py` (PR #313) -- an N-alkylated sulfonamide
  is out of scope, since 'sulfamoyl' with N-substituents needs its own
  N-prefix handling not built here.
- Otherwise mirrors `_carboxylic_acid.py`'s acyclic-chain path: P-14.3.3's
  locant-omission rule for the -COOH suffix itself (always C1, never
  cited), P-14.3.4.2(a)'s locant-omission rule for a mononuclear parent's
  own substituent locant, and 'sulfamoyl' is injected into the same
  {atom_idx -> name} map `_carboxylic_acid.py` already uses for halogens/
  a demoted hydroxyl, so it is formatted, alphabetized, and multiplied by
  the same general substituent-prefix machinery.

Scope, deliberately narrow (mirrors `_carboxylic_acid_sulfinic_acid.py`):
a single carboxylic acid plus one or more *unsubstituted* sulfonamides,
all on one acyclic saturated chain, with halogen substituents allowed.
One narrow *aromatic*-ring exception:
`_name_phenyl_chain_carboxylic_acid_sulfonamide` names a carboxylic
acid/sulfonamide chain hanging off a single plain, unsubstituted benzene
ring (e.g. '3-phenyl-2-sulfamoylpropanoic acid', PubChem CID 70062822),
mirroring `_carboxylic_acid_sulfonic_acid.py`'s identical benzene-ring-
substituent path. Explicitly out of scope (raise `UnsupportedStructure`):
any chain unsaturation (ene/yne), any ring other than the single benzene-
substituent exception above, more than one carboxylic acid, an
N-substituted sulfonamide, a coexisting standalone hydroxyl/ether/other
heteroatom, a specified stereocenter, and any carboxylic acid/sulfonamide
not captured by a single longest chain.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    longest_chains,
    lowest_locant_set,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._seniority import senior_class
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, 16, *HALOGEN_PREFIXES}

assert senior_class("carboxylic_acid", "sulfonamide") == "carboxylic_acid"


def _sulfonamide_sulfur_atoms(mol):
    """Sulfur atoms shaped like an *unsubstituted* sulfonamide group
    (-SO2NH2): bonded to exactly one carbon, two double-bonded (terminal)
    oxygens, and one single-bonded nitrogen that is itself terminal (two
    hydrogens, no other substituents). Mirrors
    `_sulfonic_acid_sulfonamide.py`'s identical helper. N-alkylated
    sulfonamides are excluded -- out of scope for this narrow module."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(oxygens) != 2 or len(nitrogens) != 1:
            continue
        double_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(double_os) != 2 or any(o.GetDegree() != 1 for o in double_os):
            continue
        (nitrogen,) = nitrogens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if nitrogen.GetDegree() != 1 or nitrogen.GetTotalNumHs() != 2:
            continue
        matches.append(atom)
    return matches


def _carboxyl_carbons(mol):
    """Carbon atoms shaped like a -COOH group: one double-bonded (terminal)
    oxygen and one single-bonded, one-H hydroxyl oxygen, both otherwise
    unaccounted-for. Mirrors the core shape check `_carboxylic_acid.py`
    uses, narrowed to plain -COOH only (no standalone-hydroxyl carve-out,
    kept out of scope here)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) != 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and o.GetTotalNumHs() == 1
        ]
        if len(carbonyls) != 1 or len(hydroxyls) != 1:
            continue
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            continue
        matches.append(atom)
    return matches


def has_carboxylic_acid_sulfonamide_shape(mol) -> bool:
    return bool(_carboxyl_carbons(mol)) and bool(_sulfonamide_sulfur_atoms(mol))


def _validate_and_collect(mol, aromatic_ring_atoms=frozenset()):
    carboxyl_carbons = _carboxyl_carbons(mol)
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "exactly one carboxylic acid is required; this module only "
            "handles a single carboxylic acid combined with one or more "
            "sulfonamides"
        )
    (carboxyl_carbon,) = carboxyl_carbons
    carboxyl_oxygens = {n.GetIdx() for n in carboxyl_carbon.GetNeighbors() if n.GetAtomicNum() == 8}

    sulfonamides = _sulfonamide_sulfur_atoms(mol)
    if not sulfonamides:
        raise UnsupportedStructure(
            "no unsubstituted sulfonamide found; this module only "
            "handles a carboxylic acid combined with at least one plain "
            "-SO2NH2 sulfonamide (see _carboxylic_acid.py for a plain "
            "carboxylic acid, and _sulfonamide.py for N-substituted forms)"
        )
    sulfonamide_idxs = {s.GetIdx() for s in sulfonamides}
    sulfonamide_oxygens = {
        o.GetIdx() for s in sulfonamides for o in s.GetNeighbors() if o.GetAtomicNum() == 8
    }
    sulfonamide_nitrogens = {
        n.GetIdx() for s in sulfonamides for n in s.GetNeighbors() if n.GetAtomicNum() == 7
    }

    accounted_oxygen_idxs = carboxyl_oxygens | sulfonamide_oxygens

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the carboxylic acid group's own "
                "oxygens, a plain sulfonamide's -NH2 nitrogen/oxygens "
                "(P-65.1.1/P-65.3.1), and halogen substituents (P-35.2.1) "
                "are not supported yet"
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
        elif atomic_num == 7:
            if atom.GetIdx() not in sulfonamide_nitrogens:
                raise UnsupportedStructure(
                    "a nitrogen that isn't part of a plain -SO2NH2 "
                    "sulfonamide is out of scope for this module (e.g. an "
                    "N-substituted sulfonamide or a coexisting amine)"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in accounted_oxygen_idxs:
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the carboxylic acid or "
                    "sulfonamide groups is out of scope for this module "
                    "(e.g. a coexisting hydroxyl, ether, or carbonyl)"
                )
        elif atomic_num == 16:
            if atom.GetIdx() not in sulfonamide_idxs:
                raise UnsupportedStructure(
                    "a sulfur atom not shaped like a plain unsubstituted "
                    "sulfonamide is out of scope for this module"
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
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    return carboxyl_carbon.GetIdx(), carboxyl_oxygens, sulfonamide_idxs


def _name_from_substituents(chain_length, grouped):
    # P-14.3.4.2(a): a mononuclear parent's substituent locant is always
    # '1' and never cited.
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    stem = alkane_name(chain_length)
    # P-14.3.3: the -COOH suffix locant is never cited, since a -COOH
    # carbon is always the chain-terminal C1.
    return f"{prefix}{stem[:-1]}oic acid"


def _candidate_key(chain_length, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, grouped)
    return (locant_set, citation_locants, name), name


def _substituents_for_chain(graph, chain, names, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, names) for root in branch_roots]
    return substituents


def _name_phenyl_chain_carboxylic_acid_sulfonamide(mol, ring_atoms):
    """Name a carboxylic acid plus one or more unsubstituted sulfonamides,
    all lying on a single unbranched chain hanging off one atom of an
    otherwise-plain, unsubstituted benzene ring -- e.g.
    '3-phenyl-2-sulfamoylpropanoic acid' (PubChem CID 70062822). The ring
    is cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_carboxylic_acid_sulfonic_acid.py`'s
    `_name_phenyl_chain_carboxylic_acid_sulfonic_acid` (PR #362) with the
    demoted group swapped from sulfonic acid/'sulfo' to
    sulfonamide/'sulfamoyl'. Narrower than the acyclic path above: no
    chain unsaturation and no specified stereocenter."""
    carboxyl_carbon, carboxyl_oxygens, sulfonamide_idxs = _validate_and_collect(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "carboxylic acid/sulfonamide chain is not supported yet"
        )
    excluded_from_unsaturation_check = {carboxyl_carbon} | sulfonamide_idxs
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_from_unsaturation_check
        and b[1] not in excluded_from_unsaturation_check
        and b[0] not in ring_atoms
        and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "carboxylic acid/sulfonamide chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain carboxylic acid/sulfonamide is not "
            "supported yet"
        )
    ring_atom, chain_root = attachment
    chain, _ = longest_branched_chain(graph, carboxyl_carbon, ring_atoms, carboxyl_oxygens | sulfonamide_idxs)
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a -COOH group directly attached to the benzene ring (no "
            "intervening chain carbon) uses the separate 'carboxylic "
            "acid' suffix construction (P-65.1.1.2), out of scope for "
            "this acyclic-chain-parent module"
        )
    chain_set = set(chain)
    sulfonamide_carbons = {
        n.GetIdx()
        for s in sulfonamide_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    if not sulfonamide_carbons.issubset(chain_set):
        raise UnsupportedStructure(
            "a sulfonamide outside the single unbranched chain hanging "
            "off the benzene ring is not supported yet"
        )

    chain_length = len(chain)
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    names = {**halogen_substituents(mol), **{s: "sulfamoyl" for s in sulfonamide_idxs}}
    substituents = _substituents_for_chain(graph, chain, names, carboxyl_oxygens | {ring_atom})
    ring_entry = name_branch(graph, ring_atom, chain_root, names, ring_atoms)
    substituents.setdefault(position_of[chain_root], []).append(ring_entry)
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, grouped)


def name_carboxylic_acid_sulfonamide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_carboxylic_acid_sulfonamide(mol, ring_atoms)

    carboxyl_carbon, carboxyl_oxygens, sulfonamide_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a carboxylic acid/sulfonamide combination on a ring is out "
            "of scope for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    excluded_from_unsaturation_check = {carboxyl_carbon} | sulfonamide_idxs
    if any(
        a not in excluded_from_unsaturation_check and b not in excluded_from_unsaturation_check
        for a, b, _ in all_non_single
    ):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a carboxylic acid/"
            "sulfonamide combination is out of scope for this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{s: "sulfamoyl" for s in sulfonamide_idxs}}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    sulfonamide_carbons = {
        n.GetIdx()
        for s in sulfonamide_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    eligible = [
        chain
        for chain in chains
        if carboxyl_carbon in chain and sulfonamide_carbons.issubset(set(chain))
    ]
    if not eligible:
        raise UnsupportedStructure(
            "not every carboxylic acid/sulfonamide-bearing carbon lies on "
            "a single longest carbon chain; a shorter principal chain is "
            "not supported yet"
        )

    excluded = carboxyl_oxygens
    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != carboxyl_carbon:
                # A -COOH carbon must sit at C1 (P-65.1.1); a direction
                # that doesn't start there is never valid.
                continue
            substituents = _substituents_for_chain(graph, candidate, names, excluded)
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
