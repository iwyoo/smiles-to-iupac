"""Naming of a molecule combining exactly one sulfonic acid (-SO3H) with
one or more sulfinic acid (-S(=O)OH) groups on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41/P-43 (`_seniority.py`): sulfonic acids (Table 4.4 class 4) outrank
  sulfinic acids (class 9), so a coexisting sulfinic acid is demoted to
  the 'sulfino' substituent prefix (P-65.3.1) instead of its own
  '-sulfinic acid' suffix -- e.g. 'OS(=O)CCS(=O)(=O)O' ->
  '2-sulfinoethane-1-sulfonic acid'. This mirrors
  `_sulfonic_acid_thiol.py`, the first module built on `_seniority.py`,
  with the demoted group swapped from thiol/'sulfanyl' to sulfinic
  acid/'sulfino'.
- Otherwise reuses `_sulfonic_acid.py`'s own acyclic-chain naming function
  directly (`_name_acyclic_sulfonic_acid`, via
  `_coexisting_groups.name_via_senior_acyclic`) rather than duplicating
  its chain-search/numbering logic, mirroring the
  `_sulfonic_acid_thiol.py` pilot (PR #413) -- 'sulfino' is injected as
  an `extra_names` prefix, and the sulfinic-bearing carbon(s) are passed
  as `required_atoms`. The P-14.3.4.2(a)/(b) locant-omission rules and
  the -SO3H locant minimized before substituent-prefix locants (P-45.2)
  both come from that shared function unchanged.

Scope, deliberately narrow (mirrors `_sulfonic_acid_thiol.py`): a single
sulfonic acid plus one or more sulfinic acids, all on one acyclic
saturated chain, with halogen substituents allowed. Explicitly out of
scope (raise `UnsupportedStructure`): any chain unsaturation (ene/yne),
any ring, more than one sulfonic acid, a coexisting hydroxyl/ether/other
heteroatom, a specified stereocenter, and any sulfonic acid/sulfinic acid
not captured by a single longest chain.
"""

from rdkit import Chem

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    non_single_bonds,
)
from ._sulfonic_acid import _name_acyclic_sulfonic_acid

_ALLOWED_ATOMIC_NUMS = {6, 8, 16, *HALOGEN_PREFIXES}


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


def has_sulfonic_acid_sulfinic_acid_shape(mol) -> bool:
    return bool(_sulfonic_sulfur_atoms(mol)) and bool(_sulfinic_sulfur_atoms(mol))


def _validate_and_collect(mol):
    sulfonic_atoms = _sulfonic_sulfur_atoms(mol)
    if len(sulfonic_atoms) != 1:
        raise UnsupportedStructure(
            "exactly one sulfonic acid is required; this module only "
            "handles a single sulfonic acid combined with one or more "
            "sulfinic acids"
        )
    (sulfonic_sulfur,) = sulfonic_atoms
    (so3h_carbon,) = (n for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 6)
    sulfonic_oxygens = {n.GetIdx() for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 8}

    sulfinics = _sulfinic_sulfur_atoms(mol)
    if not sulfinics:
        raise UnsupportedStructure(
            "no sulfinic acid found; this module only handles a sulfonic "
            "acid combined with at least one sulfinic acid (see "
            "_sulfonic_acid.py for a plain sulfonic acid)"
        )
    sulfinic_idxs = {s.GetIdx() for s in sulfinics}
    sulfinic_oxygens = {
        o.GetIdx() for s in sulfinics for o in s.GetNeighbors() if o.GetAtomicNum() == 8
    }

    accounted_sulfur_idxs = {sulfonic_sulfur.GetIdx()} | sulfinic_idxs
    accounted_oxygen_idxs = sulfonic_oxygens | sulfinic_oxygens

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the sulfonic/sulfinic acid groups' "
                "own oxygens (P-65.3.1) and halogen substituents "
                "(P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in accounted_oxygen_idxs:
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the sulfonic/sulfinic "
                    "acid groups is out of scope for this module (e.g. a "
                    "coexisting hydroxyl, ether, or carbonyl)"
                )
        elif atomic_num == 16:
            if atom.GetIdx() not in accounted_sulfur_idxs:
                raise UnsupportedStructure(
                    "a sulfur atom not shaped like the sulfonic acid group "
                    "or a plain sulfinic acid is out of scope for this module"
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

    return sulfonic_sulfur.GetIdx(), so3h_carbon.GetIdx(), sulfinic_idxs


def name_sulfonic_acid_sulfinic_acid(mol) -> str:
    sulfonic_sulfur_idx, so3h_carbon, sulfinic_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a sulfonic acid/sulfinic acid combination on a ring is out "
            "of scope for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    acid_sulfur_idxs = {sulfonic_sulfur_idx} | sulfinic_idxs
    if any(a not in acid_sulfur_idxs and b not in acid_sulfur_idxs for a, b, _ in all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a sulfonic acid/"
            "sulfinic acid combination is out of scope for this module"
        )

    sulfinic_carbons = {
        n.GetIdx()
        for s in sulfinic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    return name_via_senior_acyclic(
        _name_acyclic_sulfonic_acid,
        "sulfonic_acid",
        "sulfinic_acid",
        (mol, sulfonic_sulfur_idx, so3h_carbon, []),
        {s: "sulfino" for s in sulfinic_idxs},
        required_atoms=sulfinic_carbons,
    )
