"""Naming of a molecule combining one or more thiols (-SH) with exactly one
primary amine (-NH2) on the same acyclic saturated carbon chain, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 (`_seniority.py`): '-thiol' (`_thiol.py`, Table 4.4 class 49 --
  shared with alcohol as chalcogen analogues) far outranks 'amine' (class
  51), so a coexisting primary amine is demoted to the 'amino' substituent
  prefix instead of its own '-amine' suffix, e.g. 'NCCS' (cysteamine) ->
  '2-aminoethane-1-thiol' (PubChem CID 5726's own IUPACName is
  '2-aminoethanethiol', the retained 'ethanethiol' stem kept even when
  substituted; this module follows `_alcohol_amine.py`'s already-
  established divergence from that convention instead, mirroring
  `_thiol.py`'s own '-1-thiol' locant citation once any substituent is
  present).
  Reuses `_thiol.py`'s own `_name_acyclic_thiol` chain-search/numbering
  function directly (via `_coexisting_groups.py`), the same pattern
  `_alcohol_amine.py` used for `_alcohol.py`.
- Otherwise mirrors `_thiol.py`'s acyclic path exactly: the -SH locant(s)
  are chosen to be as low as possible (P-44), and the amine nitrogen's own
  locant is always cited via the 'amino' prefix.

Scope, deliberately narrow (mirrors `_alcohol_amine.py`): one or more
thiols plus a single primary amine, both on one acyclic *saturated*
chain, with halogen substituents allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any specified stereocenter, any ring, more than
one amine, a secondary/tertiary amine, a coexisting heteroatom other than
the thiols/amine/halogens, and any thiol/amine not captured by a single
longest chain. A benzene-ring-substituent chain variant (mirroring
`_thiol.py`'s `_name_phenyl_chain_thiol`) is a separate follow-up.
"""

from rdkit import Chem

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    adjacency,
    non_single_bonds,
    specified_stereocenters,
    validate_allowed_atoms,
)
from ._thiol import _name_acyclic_thiol



def _find_thiols(mol):
    thiols = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0 or atom.GetTotalNumHs() != 1:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            thiols.add(atom.GetIdx())
    return thiols


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


def has_thiol_amine_shape(mol) -> bool:
    return bool(_find_thiols(mol)) and len(_find_primary_amines(mol)) == 1


def _validate(mol, thiols, amines):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than a thiol sulfur (P-41), a primary amine "
        "nitrogen (P-41, Table 3.3), and halogen substituents (P-35.2.1) "
        "are not supported yet",
        [
            (
                16,
                thiols,
                "a sulfur atom that isn't a plain thiol (-SH) is out of "
                "scope for this module (e.g. a coexisting sulfide or "
                "sulfonic acid)",
            ),
            (
                7,
                amines,
                "a nitrogen that isn't a plain primary amine (-NH2) is out "
                "of scope for this module (secondary/tertiary amines, "
                "imines, and nitriles are not supported)",
            ),
        ],
    )


def name_thiol_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a thiol/amine combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    thiols = _find_thiols(mol)
    amines = _find_primary_amines(mol)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with one or more "
            "thiols is supported here (see _thiol.py for a plain thiol, "
            "_amine.py for a plain amine)"
        )
    _validate(mol, thiols, amines)

    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a thiol/amine "
            "seniority demotion is out of scope for this module (1st-pass "
            "scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a thiol/amine seniority "
            "demotion is not supported yet"
        )

    graph = adjacency(mol)
    (amine_nitrogen,) = amines
    (amine_carbon,) = graph[amine_nitrogen]
    return name_via_senior_acyclic(
        _name_acyclic_thiol,
        "alcohol",
        "amine",
        (mol, thiols, []),
        {amine_nitrogen: "amino"},
        required_atoms={amine_carbon},
    )
