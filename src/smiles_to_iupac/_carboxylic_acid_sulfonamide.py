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
- Both the acyclic-chain path and the benzene-ring-substituent path reuse
  `_carboxylic_acid.py`'s own naming functions directly
  (`_name_acyclic_carboxylic_acid`/`_name_phenyl_chain_carboxylic_acid`,
  via `_coexisting_groups.name_via_senior_acyclic`/
  `name_via_senior_phenyl_chain`) rather than duplicating their
  chain-search/numbering logic, mirroring the
  `_carboxylic_acid_sulfinic_acid.py` acyclic migration (PR #510) and the
  `_carboxylic_acid_sulfinic_acid.py`/`_carboxylic_acid_sulfonic_acid.py`
  phenyl-chain migrations (PR #513/#514) -- 'sulfamoyl' is injected as an
  `extra_names` prefix, and the sulfonamide-bearing carbon(s) are passed
  as `required_atoms`.

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

from ._coexisting_groups import name_via_senior_acyclic, name_via_senior_phenyl_chain
from ._common import (
    UnsupportedStructure,
    is_plain_benzene_ring,
    non_single_bonds,
    unsubstituted_sulfonamide_sulfur_atoms,
    validate_allowed_atoms,
)
from ._carboxylic_acid import _name_acyclic_carboxylic_acid, _name_phenyl_chain_carboxylic_acid


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
    return bool(_carboxyl_carbons(mol)) and bool(unsubstituted_sulfonamide_sulfur_atoms(mol))


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

    sulfonamides = unsubstituted_sulfonamide_sulfur_atoms(mol)
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

    validate_allowed_atoms(
        mol,
        "heteroatoms other than the carboxylic acid group's own oxygens, "
        "a plain sulfonamide's -NH2 nitrogen/oxygens (P-65.1.1/P-65.3.1), "
        "and halogen substituents (P-35.2.1) are not supported yet",
        [
            (
                7,
                sulfonamide_nitrogens,
                "a nitrogen that isn't part of a plain -SO2NH2 sulfonamide "
                "is out of scope for this module (e.g. an N-substituted "
                "sulfonamide or a coexisting amine)",
            ),
            (
                8,
                accounted_oxygen_idxs,
                "an oxygen that isn't part of the carboxylic acid or "
                "sulfonamide groups is out of scope for this module (e.g. "
                "a coexisting hydroxyl, ether, or carbonyl)",
            ),
            (
                16,
                sulfonamide_idxs,
                "a sulfur atom not shaped like a plain unsubstituted "
                "sulfonamide is out of scope for this module",
            ),
        ],
        aromatic_ring_atoms=aromatic_ring_atoms,
    )

    return carboxyl_carbon.GetIdx(), carboxyl_oxygens, sulfonamide_idxs


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
    chain unsaturation and no specified stereocenter.

    The chain-search/numbering itself reuses `_carboxylic_acid.py`'s own
    `_name_phenyl_chain_carboxylic_acid` (via `_coexisting_groups.
    name_via_senior_phenyl_chain`), same as the acyclic path above and
    mirroring the `_carboxylic_acid_sulfinic_acid.py`/
    `_carboxylic_acid_sulfonic_acid.py` migrations (PR #513/#514) --
    'sulfamoyl' is injected as an `extra_names` prefix, and the
    sulfonamide-bearing carbon(s)/oxygens/nitrogen are passed as
    `required_atoms`/`extra_accounted_atoms` (the nitrogen must be
    included here too, unlike the sulfinic/sulfonic acid cases, since a
    plain sulfonamide group has one)."""
    _carboxyl_carbon, _carboxyl_oxygens, sulfonamide_idxs = _validate_and_collect(
        mol, aromatic_ring_atoms=ring_atoms
    )

    sulfonamide_carbons = {
        n.GetIdx()
        for s in sulfonamide_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    sulfonamide_other_atoms = {
        n.GetIdx()
        for s in sulfonamide_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() in (7, 8)
    }
    return name_via_senior_phenyl_chain(
        _name_phenyl_chain_carboxylic_acid,
        "carboxylic_acid",
        "sulfonamide",
        mol,
        ring_atoms,
        {s: "sulfamoyl" for s in sulfonamide_idxs},
        required_atoms=sulfonamide_carbons,
        extra_accounted_atoms=sulfonamide_idxs | sulfonamide_other_atoms,
    )


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

    sulfonamide_carbons = {
        n.GetIdx()
        for s in sulfonamide_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    return name_via_senior_acyclic(
        _name_acyclic_carboxylic_acid,
        "carboxylic_acid",
        "sulfonamide",
        (mol, {carboxyl_carbon}, carboxyl_oxygens, set(), []),
        {s: "sulfamoyl" for s in sulfonamide_idxs},
        required_atoms=sulfonamide_carbons,
    )
