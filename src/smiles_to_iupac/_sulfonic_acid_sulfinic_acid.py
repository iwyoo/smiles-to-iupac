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

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    non_single_bonds,
    validate_allowed_atoms,
)
from ._sulfinic_acid import _sulfinic_sulfur_atoms
from ._sulfonic_acid import _name_acyclic_sulfonic_acid, _sulfonic_sulfur_atoms


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

    validate_allowed_atoms(
        mol,
        "heteroatoms other than the sulfonic/sulfinic acid groups' own "
        "oxygens (P-65.3.1) and halogen substituents (P-35.2.1) are not "
        "supported yet",
        [
            (
                8,
                accounted_oxygen_idxs,
                "an oxygen that isn't part of the sulfonic/sulfinic acid "
                "groups is out of scope for this module (e.g. a "
                "coexisting hydroxyl, ether, or carbonyl)",
            ),
            (
                16,
                accounted_sulfur_idxs,
                "a sulfur atom not shaped like the sulfonic acid group or "
                "a plain sulfinic acid is out of scope for this module",
            ),
        ],
    )

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
