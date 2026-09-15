"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
one or more sulfinic acid (-S(=O)OH) groups on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41/P-43 (`_seniority.py`): carboxylic acid (Table 4.4 class 1)
  outranks sulfinic acid (class 9), so a coexisting sulfinic acid is
  demoted to the 'sulfino' substituent prefix (P-65.3.1) instead of its
  own '-sulfinic acid' suffix -- e.g. 'OC(=O)CS(=O)O' -> '2-sulfinoacetic
  acid' (PubChem CID 13101808 confirms the systematic '2-sulfinoacetic
  acid' form; this module always uses the systematic '...anoic acid'
  stem rather than a retained name like 'acetic acid'/'formic acid',
  mirroring `_carboxylic_acid.py`'s and
  `_carboxylic_acid_sulfonic_acid.py`'s own established convention once
  any substituent is present). Mirrors `_carboxylic_acid_sulfonic_acid.py`
  (PR #314), the first pairwise module where carboxylic acid wins the
  suffix, with the demoted group swapped from sulfonic acid/'sulfo' to
  sulfinic acid/'sulfino'.
- Both the acyclic-chain path and the benzene-ring-substituent path reuse
  `_carboxylic_acid.py`'s own naming functions directly
  (`_name_acyclic_carboxylic_acid`/`_name_phenyl_chain_carboxylic_acid`,
  via `_coexisting_groups.name_via_senior_acyclic`/
  `name_via_senior_phenyl_chain`) rather than duplicating their
  chain-search/numbering logic, mirroring the
  `_carboxylic_acid_seleninic_acid.py` acyclic migration (PR #507) --
  'sulfino' is injected as an `extra_names` prefix, and the
  sulfinic-bearing carbon(s) are passed as `required_atoms`.

Scope, deliberately narrow (mirrors `_carboxylic_acid_sulfonic_acid.py`):
a single carboxylic acid plus one or more sulfinic acids, all on one
acyclic saturated chain, with halogen substituents allowed. One narrow
*aromatic*-ring exception: `_name_phenyl_chain_carboxylic_acid_sulfinic_acid`
names a carboxylic acid/sulfinic acid chain hanging off a single plain,
unsubstituted benzene ring (e.g. '3-phenyl-3-sulfinopropanoic acid',
PubChem CID 174325886), mirroring
`_carboxylic_acid_sulfonic_acid.py`'s identical benzene-ring-substituent
path. Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring other than the single benzene-
substituent exception above, more than one carboxylic acid, a coexisting
standalone hydroxyl/ether/other heteroatom, a specified stereocenter, and
any carboxylic acid/sulfinic acid not captured by a single longest chain.
"""

from ._coexisting_groups import name_via_senior_acyclic, name_via_senior_phenyl_chain
from ._common import (
    UnsupportedStructure,
    halogen_substituents,
    is_plain_benzene_ring,
    non_single_bonds,
    validate_allowed_atoms,
)
from ._carboxylic_acid import _name_acyclic_carboxylic_acid, _name_phenyl_chain_carboxylic_acid


def _sulfinic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfinic acid group: bonded to exactly
    one carbon, one double-bonded (terminal) oxygen, and one
    single-bonded hydroxyl oxygen (terminal, one H). Mirrors
    `_sulfinic_acid.py`'s identical helper."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 2:
            continue
        double_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
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


def has_carboxylic_acid_sulfinic_acid_shape(mol) -> bool:
    return bool(_carboxyl_carbons(mol)) and bool(_sulfinic_sulfur_atoms(mol))


def _validate_and_collect(mol, aromatic_ring_atoms=frozenset()):
    carboxyl_carbons = _carboxyl_carbons(mol)
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "exactly one carboxylic acid is required; this module only "
            "handles a single carboxylic acid combined with one or more "
            "sulfinic acids"
        )
    (carboxyl_carbon,) = carboxyl_carbons
    carboxyl_oxygens = {n.GetIdx() for n in carboxyl_carbon.GetNeighbors() if n.GetAtomicNum() == 8}

    sulfinics = _sulfinic_sulfur_atoms(mol)
    if not sulfinics:
        raise UnsupportedStructure(
            "no sulfinic acid found; this module only handles a "
            "carboxylic acid combined with at least one sulfinic acid "
            "(see _carboxylic_acid.py for a plain carboxylic acid)"
        )
    sulfinic_idxs = {s.GetIdx() for s in sulfinics}
    sulfinic_oxygens = {
        o.GetIdx() for s in sulfinics for o in s.GetNeighbors() if o.GetAtomicNum() == 8
    }

    accounted_oxygen_idxs = carboxyl_oxygens | sulfinic_oxygens

    validate_allowed_atoms(
        mol,
        "heteroatoms other than the carboxylic/sulfinic acid groups' own "
        "oxygens (P-65.1.1/P-65.3.1) and halogen substituents (P-35.2.1) "
        "are not supported yet",
        [
            (
                8,
                accounted_oxygen_idxs,
                "an oxygen that isn't part of the carboxylic/sulfinic acid "
                "groups is out of scope for this module (e.g. a "
                "coexisting hydroxyl, ether, or carbonyl)",
            ),
            (
                16,
                sulfinic_idxs,
                "a sulfur atom not shaped like a sulfinic acid group is "
                "out of scope for this module",
            ),
        ],
        aromatic_ring_atoms=aromatic_ring_atoms,
    )

    return carboxyl_carbon.GetIdx(), carboxyl_oxygens, sulfinic_idxs


def _name_phenyl_chain_carboxylic_acid_sulfinic_acid(mol, ring_atoms):
    """Name a carboxylic acid plus one or more sulfinic acids, all lying
    on a single unbranched chain hanging off one atom of an otherwise-
    plain, unsubstituted benzene ring -- e.g. '3-phenyl-3-sulfinopropanoic
    acid' (PubChem CID 174325886). The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring
    `_carboxylic_acid_sulfonic_acid.py`'s
    `_name_phenyl_chain_carboxylic_acid_sulfonic_acid` (PR #362) with the
    demoted group swapped from sulfonic acid/'sulfo' to sulfinic
    acid/'sulfino'. Narrower than the acyclic path above: no chain
    unsaturation and no specified stereocenter.

    The chain-search/numbering itself reuses `_carboxylic_acid.py`'s own
    `_name_phenyl_chain_carboxylic_acid` (via `_coexisting_groups.
    name_via_senior_phenyl_chain`), same as the acyclic path above --
    'sulfino' is injected as an `extra_names` prefix, and the
    sulfinic-bearing carbon(s)/oxygens are passed as `required_atoms`/
    `extra_accounted_atoms`. The specified-stereocenter and chain-
    unsaturation rejections are the shared function's own (same checks,
    now unduplicated)."""
    _carboxyl_carbon, _carboxyl_oxygens, sulfinic_idxs = _validate_and_collect(
        mol, aromatic_ring_atoms=ring_atoms
    )

    sulfinic_carbons = {
        n.GetIdx()
        for s in sulfinic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    sulfinic_oxygens = {
        n.GetIdx()
        for s in sulfinic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 8
    }
    return name_via_senior_phenyl_chain(
        _name_phenyl_chain_carboxylic_acid,
        "carboxylic_acid",
        "sulfinic_acid",
        mol,
        ring_atoms,
        {s: "sulfino" for s in sulfinic_idxs},
        required_atoms=sulfinic_carbons,
        extra_accounted_atoms=sulfinic_idxs | sulfinic_oxygens,
    )


def name_carboxylic_acid_sulfinic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_carboxylic_acid_sulfinic_acid(mol, ring_atoms)

    carboxyl_carbon, carboxyl_oxygens, sulfinic_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a carboxylic acid/sulfinic acid combination on a ring is out "
            "of scope for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    excluded_from_unsaturation_check = {carboxyl_carbon} | sulfinic_idxs
    if any(
        a not in excluded_from_unsaturation_check and b not in excluded_from_unsaturation_check
        for a, b, _ in all_non_single
    ):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a carboxylic acid/"
            "sulfinic acid combination is out of scope for this module"
        )

    sulfinic_carbons = {
        n.GetIdx()
        for s in sulfinic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    return name_via_senior_acyclic(
        _name_acyclic_carboxylic_acid,
        "carboxylic_acid",
        "sulfinic_acid",
        (mol, {carboxyl_carbon}, carboxyl_oxygens, set(), []),
        {s: "sulfino" for s in sulfinic_idxs},
        required_atoms=sulfinic_carbons,
    )
