"""Naming of thiocyanates (R-S-C#N), per the IUPAC 2013 Recommendations
("the Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a thiocyanate
  ester is named the same way any other retained-acid ester is: 'R
  thiocyanate', two words, the alkyl group R cited first exactly as
  `_carbamate.py`'s/`_ester.py`'s alcohol part is, e.g. 'ethyl thiocyanate'
  for CH3CH2-S-C#N. The acid side is fixed: it is always the single word
  'thiocyanate', with no chain-length logic of its own -- mirroring
  `_carbamate.py`'s own 'carbamate' word exactly.
- This is the constitutional isomer of isothiocyanate (R-N=C=S,
  `_isothiocyanate.py`, already covered): the same R-S-C#N vs. R-N=C=S
  distinction the Blue Book itself draws between the two classes.
- Confirmed via PubChem structure match: `CSC#N` -> "methyl thiocyanate",
  `CCSC#N` -> "ethyl thiocyanate", `CCCSC#N` -> "propyl thiocyanate". The
  chalcogen analogues (cyanate, selenocyanate, tellurocyanate) also exist
  as their own PubChem-confirmed structures but are deferred to a
  follow-up task, mirroring how the isocyanate/isothiocyanate/
  isoselenocyanate/isotellurocyanate family was built up one module at a
  time.

Scope, deliberately narrow (mirrors `_carbamate.py`'s own first pass): R
is restricted to a plain, unsubstituted, saturated, acyclic alkyl group
(branched or unbranched) attached at its own chain terminus (e.g.
'methyl', 'propan-2-yl', 'tert-butyl') -- R's name is built with
`name_branch` (P-29 PIN style, fixed project-wide by PR #237; mirrors
`_carbamate.py`'s/`_ester.py`'s identical fix, PR #328/#329), never
parenthesized (the "R thiocyanate" two-word pattern has no nested-prefix
ambiguity to guard against). Confirmed via PubChem: `CC(C)SC#N` ->
"propan-2-yl thiocyanate" (CID 246911), `CC(C)(C)SC#N` -> "tert-butyl
thiocyanate" (CID 641640). R may also be a plain, unsubstituted benzene
ring bonded directly to the thiocyanate sulfur, e.g. 'phenyl thiocyanate'
(CID 21357) -- a chain spacer between the ring and the sulfur (e.g.
'benzyl thiocyanate') is still deferred, mirroring `_cyanate.py`'s
identical scope note. A substituted, unsaturated, or otherwise
ring-bearing R is still deferred.
Explicitly out of scope (raise `UnsupportedStructure`): any ring other
than the single plain-benzene-bonded-directly-to-S exception above, more
than one thiocyanate group, and any other heteroatom/oxygen not part of
this single thiocyanate group.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, is_plain_benzene_ring, non_single_bonds
from ._substituents import name_branch

_YNE_ORDER = 3.0


def _thiocyanate_cores(mol):
    """List of (sulfur_idx, alkyl_carbon_idx) for every R-S-C#N pattern: a
    sulfur bonded to exactly one carbon (R) by a single bond and one
    carbon (the nitrile carbon) by a single bond, that nitrile carbon
    triple-bonded to a terminal nitrogen and bonded to nothing else."""
    cores = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 2:
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


def has_thiocyanate_shape(mol) -> bool:
    return bool(_thiocyanate_cores(mol))


def name_thiocyanate(mol) -> str:
    cores = _thiocyanate_cores(mol)
    if len(cores) != 1:
        raise UnsupportedStructure(
            "exactly one thiocyanate (-S-C#N) group is required; zero or "
            "multiple such groups are not supported yet"
        )
    sulfur_idx, nitrile_c_idx, alkyl_c_idx = cores[0]

    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    if ring_info.NumRings() == 1:
        candidate_ring_atoms = set(ring_info.AtomRings()[0])
        if not is_plain_benzene_ring(mol, candidate_ring_atoms):
            raise UnsupportedStructure(
                "a non-benzene ring is out of scope for this module"
            )
        if alkyl_c_idx not in candidate_ring_atoms:
            raise UnsupportedStructure(
                "a benzene ring reached through a chain spacer (rather "
                "than bonded directly to the thiocyanate sulfur) is not "
                "supported yet"
            )
        ring_atoms = candidate_ring_atoms
    elif ring_info.NumRings() > 1:
        raise UnsupportedStructure(
            "more than one ring is out of scope for this module"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (nitrogen_idx,) = (
        n.GetIdx() for n in mol.GetAtomWithIdx(nitrile_c_idx).GetNeighbors() if n.GetIdx() != sulfur_idx
    )
    excluded = {sulfur_idx, nitrile_c_idx, nitrogen_idx}
    for atom in mol.GetAtoms():
        if atom.GetIdx() in excluded:
            continue
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms other than this single thiocyanate group are "
                "not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetIsAromatic() and atom.GetIdx() not in ring_atoms:
            raise UnsupportedStructure("aromatic rings are out of scope for this module")

    non_single = [
        b
        for b in non_single_bonds(mol)
        if nitrile_c_idx not in (b[0], b[1]) and (b[0] not in ring_atoms or b[1] not in ring_atoms)
    ]
    if non_single:
        raise UnsupportedStructure("unsaturation in the R group is not supported yet")

    r_name, _ = name_branch(adjacency(mol), alkyl_c_idx, sulfur_idx, {}, ring_atoms, mol=mol)
    return f"{r_name} thiocyanate"
