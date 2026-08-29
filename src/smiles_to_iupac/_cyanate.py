"""Naming of cyanates (R-O-C#N), per the IUPAC 2013 Recommendations ("the
Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a cyanate
  ester is named the same way any other retained-acid ester is: 'R
  cyanate', two words, the alkyl group R cited first exactly as
  `_carbamate.py`'s/`_thiocyanate.py`'s alcohol part is, e.g. 'ethyl
  cyanate' for CH3CH2-O-C#N. The acid side is fixed: it is always the
  single word 'cyanate', with no chain-length logic of its own --
  mirroring `_thiocyanate.py`'s own 'thiocyanate' word exactly, just with
  oxygen instead of sulfur.
- Confirmed via PubChem structure match: `COC#N` -> "methyl cyanate",
  `CCOC#N` -> "ethyl cyanate", `CCCOC#N` -> "propyl cyanate". The other
  chalcogen analogues (selenocyanate, tellurocyanate) also exist as their
  own PubChem-confirmed structures but are deferred to a follow-up task,
  mirroring how the isocyanate/isothiocyanate/isoselenocyanate/
  isotellurocyanate family was built up one module at a time.

Scope, deliberately narrow (mirrors `_thiocyanate.py`'s own first pass): R
is restricted to a plain, unbranched, unsubstituted, saturated alkyl group
attached at its own chain terminus (e.g. 'methyl', 'ethyl', 'propyl'); a
branched, substituted, unsaturated, or ring-bearing R is deferred.
Explicitly out of scope (raise `UnsupportedStructure`): any ring anywhere
in the molecule, more than one cyanate group, and any other
heteroatom/oxygen not part of this single cyanate group.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, carbon_adjacency, linear_branch, non_single_bonds
from ._numerals import alkyl_name

_YNE_ORDER = 3.0


def _cyanate_cores(mol):
    """List of (oxygen_idx, nitrile_carbon_idx, alkyl_carbon_idx) for every
    R-O-C#N pattern: an oxygen bonded to exactly one carbon (R) by a single
    bond and one carbon (the nitrile carbon) by a single bond, that
    nitrile carbon triple-bonded to a terminal nitrogen and bonded to
    nothing else."""
    cores = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 2:
            continue
        neighbors = atom.GetNeighbors()
        if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in neighbors):
            continue
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        if len(carbons) != 2:
            continue
        for nitrile_c, alkyl_c in ((carbons[0], carbons[1]), (carbons[1], carbons[0])):
            if nitrile_c.GetDegree() != 2:
                continue
            other_neighbors = [n for n in nitrile_c.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
            if len(other_neighbors) != 1:
                continue
            (nitrogen,) = other_neighbors
            if nitrogen.GetAtomicNum() != 7 or nitrogen.GetDegree() != 1:
                continue
            if mol.GetBondBetweenAtoms(nitrile_c.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != _YNE_ORDER:
                continue
            if nitrogen.GetFormalCharge() != 0 or nitrogen.GetIsotope() != 0:
                continue
            cores.append((atom.GetIdx(), nitrile_c.GetIdx(), alkyl_c.GetIdx()))
    return cores


def has_cyanate_shape(mol) -> bool:
    return bool(_cyanate_cores(mol))


def name_cyanate(mol) -> str:
    cores = _cyanate_cores(mol)
    if len(cores) != 1:
        raise UnsupportedStructure(
            "exactly one cyanate (-O-C#N) group is required; zero or "
            "multiple such groups are not supported yet"
        )
    oxygen_idx, nitrile_c_idx, alkyl_c_idx = cores[0]

    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-attached cyanate is out of scope for this acyclic-only "
            "module"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (nitrogen_idx,) = (
        n.GetIdx() for n in mol.GetAtomWithIdx(nitrile_c_idx).GetNeighbors() if n.GetIdx() != oxygen_idx
    )
    excluded = {oxygen_idx, nitrile_c_idx, nitrogen_idx}
    for atom in mol.GetAtoms():
        if atom.GetIdx() in excluded:
            continue
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms other than this single cyanate group are not "
                "supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0 or atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic rings are out of scope for this module")

    non_single = [b for b in non_single_bonds(mol) if nitrile_c_idx not in (b[0], b[1])]
    if non_single:
        raise UnsupportedStructure("unsaturation in the R group is not supported yet")

    carbon_graph = carbon_adjacency(mol)
    length = linear_branch(carbon_graph, alkyl_c_idx, None)
    if length is None:
        raise UnsupportedStructure("a branched R group is not supported yet")

    return f"{alkyl_name(length)} cyanate"
