"""Naming of a molecule combining exactly one sulfonic acid (-SO3H) with
one or more unsubstituted sulfonamide (-SO2NH2) groups on the same
acyclic saturated carbon chain, per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-41/P-43 (`_seniority.py`): sulfonic acids (Table 4.4 class 4) outrank
  sulfonamides (class 19), so a coexisting sulfonamide is demoted to the
  'sulfamoyl' substituent prefix (P-65.3.1) instead of its own
  '-sulfonamide' suffix -- e.g. 'NS(=O)(=O)CCS(=O)(=O)O' ->
  '2-sulfamoylethanesulfonic acid'. This mirrors
  `_sulfonic_acid_sulfinic_acid.py`, with the demoted group swapped from
  sulfinic acid/'sulfino' to sulfonamide/'sulfamoyl'.
- Otherwise reuses `_sulfonic_acid.py`'s own acyclic-chain naming function
  directly (`_name_acyclic_sulfonic_acid`, via
  `_coexisting_groups.name_via_senior_acyclic`) rather than duplicating
  its chain-search/numbering logic, mirroring the
  `_sulfonic_acid_sulfinic_acid.py` migration (PR #508) -- 'sulfamoyl' is
  injected as an `extra_names` prefix, and the sulfonamide-bearing
  carbon(s) are passed as `required_atoms`. The P-14.3.4.2(a)/(b)
  locant-omission rules and the -SO3H locant minimized before
  substituent-prefix locants (P-45.2) both come from that shared function
  unchanged.

Scope, deliberately narrow (mirrors `_sulfonic_acid_sulfinic_acid.py`): a
single sulfonic acid plus one or more *unsubstituted* sulfonamides
(-SO2NH2 only -- an N-alkylated sulfonamide is out of scope, since
'sulfamoyl' with N-substituents needs its own N-prefix handling not
built here), all on one acyclic saturated chain, with halogen
substituents allowed. Explicitly out of scope (raise
`UnsupportedStructure`): any chain unsaturation (ene/yne), any ring, more
than one sulfonic acid, an N-substituted sulfonamide, a coexisting
hydroxyl/ether/other heteroatom, a specified stereocenter, and any
sulfonic acid/sulfonamide not captured by a single longest chain.
"""

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    non_single_bonds,
    validate_allowed_atoms,
)
from ._sulfonic_acid import _name_acyclic_sulfonic_acid


def _sulfonic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfonic acid group: bonded to exactly
    one carbon, two double-bonded (terminal) oxygens, and one
    single-bonded hydroxyl oxygen (terminal, one H). Mirrors
    `_sulfonic_acid.py`'s identical helper."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 3:
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
        if len(double_os) != 2 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def _sulfonamide_sulfur_atoms(mol):
    """Sulfur atoms shaped like an *unsubstituted* sulfonamide group
    (-SO2NH2): bonded to exactly one carbon, two double-bonded (terminal)
    oxygens, and one single-bonded nitrogen that is itself terminal (two
    hydrogens, no other substituents). N-alkylated sulfonamides are
    excluded -- out of scope for this narrow module."""
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


def has_sulfonic_acid_sulfonamide_shape(mol) -> bool:
    return bool(_sulfonic_sulfur_atoms(mol)) and bool(_sulfonamide_sulfur_atoms(mol))


def _validate_and_collect(mol):
    sulfonic_atoms = _sulfonic_sulfur_atoms(mol)
    if len(sulfonic_atoms) != 1:
        raise UnsupportedStructure(
            "exactly one sulfonic acid is required; this module only "
            "handles a single sulfonic acid combined with one or more "
            "sulfonamides"
        )
    (sulfonic_sulfur,) = sulfonic_atoms
    (so3h_carbon,) = (n for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 6)
    sulfonic_oxygens = {n.GetIdx() for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 8}

    sulfonamides = _sulfonamide_sulfur_atoms(mol)
    if not sulfonamides:
        raise UnsupportedStructure(
            "no unsubstituted sulfonamide found; this module only "
            "handles a sulfonic acid combined with at least one plain "
            "-SO2NH2 sulfonamide (see _sulfonic_acid.py for a plain "
            "sulfonic acid, and _sulfonamide.py for N-substituted forms)"
        )
    sulfonamide_idxs = {s.GetIdx() for s in sulfonamides}
    sulfonamide_oxygens = {
        o.GetIdx() for s in sulfonamides for o in s.GetNeighbors() if o.GetAtomicNum() == 8
    }
    sulfonamide_nitrogens = {
        n.GetIdx() for s in sulfonamides for n in s.GetNeighbors() if n.GetAtomicNum() == 7
    }

    accounted_sulfur_idxs = {sulfonic_sulfur.GetIdx()} | sulfonamide_idxs
    accounted_oxygen_idxs = sulfonic_oxygens | sulfonamide_oxygens

    validate_allowed_atoms(
        mol,
        "heteroatoms other than the sulfonic acid group's own oxygens, a "
        "plain sulfonamide's -NH2 nitrogen/oxygens (P-65.3.1), and "
        "halogen substituents (P-35.2.1) are not supported yet",
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
                "an oxygen that isn't part of the sulfonic acid or "
                "sulfonamide groups is out of scope for this module (e.g. "
                "a coexisting hydroxyl, ether, or carbonyl)",
            ),
            (
                16,
                accounted_sulfur_idxs,
                "a sulfur atom not shaped like the sulfonic acid group or "
                "a plain unsubstituted sulfonamide is out of scope for "
                "this module",
            ),
        ],
    )

    return sulfonic_sulfur.GetIdx(), so3h_carbon.GetIdx(), sulfonamide_idxs


def name_sulfonic_acid_sulfonamide(mol) -> str:
    sulfonic_sulfur_idx, so3h_carbon, sulfonamide_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a sulfonic acid/sulfonamide combination on a ring is out of "
            "scope for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    acid_sulfur_idxs = {sulfonic_sulfur_idx} | sulfonamide_idxs
    if any(a not in acid_sulfur_idxs and b not in acid_sulfur_idxs for a, b, _ in all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a sulfonic acid/"
            "sulfonamide combination is out of scope for this module"
        )

    sulfonamide_carbons = {
        n.GetIdx()
        for s in sulfonamide_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    return name_via_senior_acyclic(
        _name_acyclic_sulfonic_acid,
        "sulfonic_acid",
        "sulfonamide",
        (mol, sulfonic_sulfur_idx, so3h_carbon, []),
        {s: "sulfamoyl" for s in sulfonamide_idxs},
        required_atoms=sulfonamide_carbons,
    )
