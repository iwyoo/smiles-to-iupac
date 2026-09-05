"""Naming of a molecule combining one or more hydroxyls (-OH) with exactly
one primary amine (-NH2) on the same acyclic saturated carbon chain, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 (`_seniority.py`): '-ol' (`_alcohol.py`, Table 4.4 class 49) far
  outranks 'amine' (class 51), so a coexisting primary amine is demoted to
  the 'amino' substituent prefix instead of its own '-amine' suffix, e.g.
  'NCCO' (ethanolamine) -> '2-aminoethan-1-ol' (PubChem CID 700's own
  IUPACName is '2-aminoethanol', the retained 'ethanol' stem kept even
  when substituted; this project's own `_alcohol.py` already diverges from
  that convention -- a substituted ethanol always cites the '-1-ol' locant,
  e.g. its own test '2-ethoxyethan-1-ol' -- so this module follows that
  same pre-existing convention rather than introduce a new divergence).
  Reuses `_alcohol.py`'s own `_name_acyclic_alcohol` chain-search/
  numbering function directly (via `_coexisting_groups.py`), injecting the
  demoted amine nitrogen as an 'amino' name the same way that function
  already injects a demoted alkoxy ether.
- Otherwise mirrors `_alcohol.py`'s acyclic path exactly: the -OH
  locant(s) are chosen to be as low as possible (P-44), and the amine
  nitrogen's own locant is always cited via the 'amino' prefix.

Scope, deliberately narrow (first pairwise pilot on `_alcohol.py`): one or
more hydroxyls plus a single primary amine, both on one acyclic
*saturated* chain, with halogen substituents allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any specified stereocenter, any ring, more than
one amine, a secondary/tertiary amine, a coexisting ether/other
heteroatom, and any hydroxyl/amine not captured by a single longest chain.
A benzene-ring-substituent chain variant (mirroring `_alcohol.py`'s
`_name_phenyl_chain_alcohol`) is a separate follow-up.
"""

from rdkit import Chem

from ._alcohol import _name_acyclic_alcohol
from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    non_single_bonds,
    specified_stereocenters,
)

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _find_hydroxyls(mol):
    hydroxyls = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0 or atom.GetTotalNumHs() != 1:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            hydroxyls.add(atom.GetIdx())
    return hydroxyls


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


def has_alcohol_amine_shape(mol) -> bool:
    return bool(_find_hydroxyls(mol)) and len(_find_primary_amines(mol)) == 1


def _validate(mol, hydroxyls, amines):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a hydroxyl oxygen (P-33.2.1), a "
                "primary amine nitrogen (P-41, Table 3.3), and halogen "
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
            if atom.GetIdx() not in hydroxyls:
                raise UnsupportedStructure(
                    "an oxygen that isn't a plain hydroxyl (-OH) is out of "
                    "scope for this module (e.g. a coexisting ether or "
                    "carbonyl)"
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


def name_alcohol_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a hydroxyl/amine combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    hydroxyls = _find_hydroxyls(mol)
    amines = _find_primary_amines(mol)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with one or more "
            "hydroxyls is supported here (see _alcohol.py for a plain "
            "alcohol, _amine.py for a plain amine)"
        )
    _validate(mol, hydroxyls, amines)

    all_non_single = non_single_bonds(mol)
    if all_non_single:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a hydroxyl/amine "
            "seniority demotion is out of scope for this module (1st-pass "
            "scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a hydroxyl/amine seniority "
            "demotion is not supported yet"
        )

    graph = adjacency(mol)
    (amine_nitrogen,) = amines
    (amine_carbon,) = graph[amine_nitrogen]
    return name_via_senior_acyclic(
        _name_acyclic_alcohol,
        "alcohol",
        "amine",
        (mol, hydroxyls, []),
        {amine_nitrogen: "amino"},
        required_atoms={amine_carbon},
    )
