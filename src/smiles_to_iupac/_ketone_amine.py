"""Naming of a molecule combining a ketone (C=O, two carbon neighbors) with
exactly one primary amine (-NH2) on the same acyclic saturated carbon chain,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 (`_seniority.py`): 'one' (`_ketone.py`, Table 4.4 class 48) far
  outranks 'amine' (class 51), so a coexisting primary amine is demoted to
  the 'amino' substituent prefix instead of its own '-amine' suffix, e.g.
  'NCC(C)=O' -> '1-aminopropan-2-one' (PubChem's own IUPACName, an exact
  match -- no retained-name divergence for this class, unlike the
  aldehyde/acid pairwise modules).
  Reuses `_ketone.py`'s own `_name_acyclic_ketone` chain-search/numbering
  function directly (via `_coexisting_groups.py`), the same pattern
  `_alcohol_amine.py`/`_thiol_amine.py` used for their own senior modules.
- Otherwise mirrors `_ketone.py`'s acyclic path exactly: the C=O locant is
  chosen to be as low as possible (P-44), and the amine nitrogen's own
  locant is always cited via the 'amino' prefix.

Scope, deliberately narrow (mirrors `_alcohol_amine.py`/`_thiol_amine.py`):
a single ketone plus a single primary amine, both on one acyclic
*saturated* chain, with halogen substituents allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any specified stereocenter, any ring, more than
one ketone or amine, a secondary/tertiary amine, a coexisting standalone
hydroxyl/other heteroatom, and any ketone/amine not captured by a single
longest chain. A benzene-ring-substituent chain variant (mirroring
`_ketone.py`'s `_name_phenyl_chain_ketone`) is a separate follow-up.
"""

from rdkit import Chem

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    non_single_bonds,
    specified_stereocenters,
)
from ._ketone import _name_acyclic_ketone

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _find_ketones(mol):
    ketones = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6 or carbon.GetIsAromatic():
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 2:
            ketones.add(atom.GetIdx())
    return ketones


def _find_primary_amines(mol):
    amines = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0 or atom.GetTotalNumHs() != 2:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            amines.add(atom.GetIdx())
    return amines


def has_ketone_amine_shape(mol) -> bool:
    return len(_find_ketones(mol)) == 1 and len(_find_primary_amines(mol)) == 1


def _validate(mol, ketones, amines):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a ketone carbonyl oxygen (P-33.4), "
                "a primary amine nitrogen (P-41, Table 3.3), and halogen "
                "substituents (P-35.2.1) are not supported yet"
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
            if atom.GetIdx() not in ketones:
                raise UnsupportedStructure(
                    "an oxygen that isn't the single ketone carbonyl is out "
                    "of scope for this module (e.g. a coexisting hydroxyl, "
                    "ether, or a second carbonyl)"
                )
        elif atomic_num == 7:
            if atom.GetIdx() not in amines:
                raise UnsupportedStructure(
                    "a nitrogen that isn't a plain primary amine (-NH2) is "
                    "out of scope for this module (secondary/tertiary "
                    "amines, imines, and nitriles are not supported)"
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


def name_ketone_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ketone/amine combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    ketones = _find_ketones(mol)
    amines = _find_primary_amines(mol)
    if len(ketones) != 1:
        raise UnsupportedStructure(
            "exactly one ketone coexisting with a primary amine is "
            "supported here (see _ketone.py for a plain ketone)"
        )
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "ketone is supported here (see _amine.py for a plain amine)"
        )
    _validate(mol, ketones, amines)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in ketones or b[1] in ketones]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a ketone/amine "
            "seniority demotion is out of scope for this module (1st-pass "
            "scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a ketone/amine seniority "
            "demotion is not supported yet"
        )

    graph = adjacency(mol)
    (amine_nitrogen,) = amines
    (amine_carbon,) = graph[amine_nitrogen]
    return name_via_senior_acyclic(
        _name_acyclic_ketone,
        "ketone",
        "amine",
        (mol, ketones, set(), []),
        {amine_nitrogen: "amino"},
        required_atoms={amine_carbon},
    )
