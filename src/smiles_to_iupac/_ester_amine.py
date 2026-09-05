"""Naming of a molecule combining an ester (R-CO-O-R') with exactly one
separate primary amine (-NH2) on the acyl chain, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-41 Table 4.1 (`_seniority.py`, `tmp/bluebook/P41.txt`): esters (class
  9) far outrank amines (class 19), so a coexisting primary amine is
  demoted to the 'amino' substituent prefix instead of its own '-amine'
  suffix, e.g. 'NCCC(=O)OC' -> 'methyl 3-aminopropanoate' (PubChem's own
  IUPACName, an exact match).
  Reuses `_ester.py`'s own `_name_acyl_part` chain-search/numbering
  function directly (via `_coexisting_groups.py`), the same pattern
  `_amide_amine.py`/`_ketone_amine.py` used for their own senior modules.
- Otherwise mirrors `_ester.py`'s plain acyclic path exactly: the acyl
  carbon is always chain C1 (P-14.3.3, no locant), the alcohol part (R')
  is a plain unsubstituted saturated acyclic alkyl group, and the amine
  nitrogen's own locant on the acyl chain is always cited via the 'amino'
  prefix.

Scope, deliberately narrow (mirrors `_amide_amine.py`/`_ketone_amine.py`):
a single ester group plus a single primary amine on the acyl chain, both
on one acyclic *saturated* chain, with halogen substituents on the acyl
chain allowed. The alcohol part (R') is restricted the same way
`_ester.py` restricts it (plain, unsubstituted, saturated, acyclic alkyl).
Explicitly out of scope (raise `UnsupportedStructure`): any amine on the
alcohol part (R') instead of the acyl chain, any ring (including an
aryl/cyclyl alcohol part), any chain unsaturation (ene/yne), any specified
stereocenter, more than one ester or amine, a secondary/tertiary amine, a
coexisting standalone hydroxyl/other heteroatom, and an amine bonded
directly to the acyl carbon itself (a different topology, out of scope).
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
from ._ester import _find_ester_group, _name_acyl_part, _name_alcohol_part

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _find_primary_amines(mol, exclude):
    amines = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() in exclude or atom.GetAtomicNum() != 7 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0 or atom.GetTotalNumHs() != 2:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            amines.add(atom.GetIdx())
    return amines


def has_ester_amine_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    try:
        acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = _find_ester_group(mol)
    except UnsupportedStructure:
        return False
    amines = _find_primary_amines(mol, {ester_oxygen.GetIdx(), carbonyl_oxygen.GetIdx()})
    return len(amines) == 1


def _validate(mol, ester_atoms, amine_nitrogen):
    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = ester_atoms
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the ester's own oxygens (P-65.6.3), "
                "a primary amine nitrogen (P-41, Table 4.1), and halogen "
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
            if atom.GetIdx() not in (carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()):
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the single ester's "
                    "carbonyl/ester pair is out of scope for this module "
                    "(e.g. a coexisting hydroxyl, ether, or second carbonyl)"
                )
        elif atomic_num == 7:
            if atom.GetIdx() != amine_nitrogen:
                raise UnsupportedStructure(
                    "a nitrogen that isn't a single plain primary amine "
                    "(-NH2) is out of scope for this module (secondary/"
                    "tertiary amines, imines, and nitriles are not "
                    "supported)"
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


def name_ester_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ester/amine combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = _find_ester_group(mol)
    excluded = {ester_oxygen.GetIdx(), carbonyl_oxygen.GetIdx()}
    amines = _find_primary_amines(mol, excluded)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "ester is supported here (see _amine.py for a plain amine)"
        )
    (amine_nitrogen,) = amines
    _validate(mol, (acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon), amine_nitrogen)

    graph = adjacency(mol)
    (amine_carbon,) = graph[amine_nitrogen]
    if amine_carbon == alcohol_carbon.GetIdx():
        raise UnsupportedStructure(
            "an amine on the alcohol part (R') rather than the acyl chain "
            "is not supported yet"
        )
    if amine_carbon == acyl_carbon.GetIdx():
        raise UnsupportedStructure(
            "an amine bonded directly to the ester's acyl carbon is not "
            "supported yet"
        )

    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    if all_non_single:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ester/amine "
            "seniority demotion is out of scope for this module (1st-pass "
            "scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ester/amine seniority "
            "demotion is not supported yet"
        )

    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())
    acyl_name = name_via_senior_acyclic(
        _name_acyl_part,
        "ester",
        "amine",
        (mol, acyl_carbon, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()),
        {amine_nitrogen: "amino"},
        required_atoms={amine_carbon},
    )
    return f"{alcohol_name} {acyl_name}"
