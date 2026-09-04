"""Naming of tellurocyanates (R-Te-C#N), per the IUPAC 2013 Recommendations
("the Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  tellurocyanate ester is named the same way any other retained-acid ester
  is: 'R tellurocyanate', two words, the alkyl group R cited first exactly
  as `_selenocyanate.py`'s/`_thiocyanate.py`'s alcohol part is, e.g.
  'methyl tellurocyanate' for CH3-Te-C#N. The acid side is fixed: it is
  always the single word 'tellurocyanate', with no chain-length logic of
  its own -- mirroring `_selenocyanate.py`'s own 'selenocyanate' word
  exactly, just with tellurium instead of selenium.
- Confirmed via PubChem structure match: `C[Te]C#N` -> "methyl
  tellurocyanate" -- the only tellurocyanate PubChem has registered at
  all (ethyl/propyl candidates both come back as CID 0, the same
  sparse-tellurium-data gap `_telluronic_acid.py`/`_tellurinic_acid.py`
  already documented). The chain-generalization mechanism itself is NOT
  independently re-verified for tellurium here -- it is inherited,
  unchanged, from the identical mechanism already confirmed by
  `_thiocyanate.py`/`_cyanate.py`/`_selenocyanate.py`, the same
  reduced-evidence bar `_telluronic_acid.py` first used for this exact
  situation.

Scope, deliberately narrow (mirrors `_selenocyanate.py`'s own first
pass): R is restricted to a plain, unsubstituted, saturated, acyclic
alkyl group (branched or unbranched) attached at its own chain terminus,
built with `name_branch` (P-29 PIN style, PR #237/#328/#329), never
parenthesized. Not independently PubChem-verified for a branched R here
(PubChem has no branched tellurocyanate registered, same sparse-tellurium
gap as above) -- inherited unchanged from the identical mechanism already
confirmed by `_thiocyanate.py`/`_cyanate.py`/`_selenocyanate.py`. A
substituted, unsaturated, or ring-bearing R is still deferred.
Explicitly out of scope (raise `UnsupportedStructure`): any
ring anywhere in the molecule, more than one tellurocyanate group, and
any other heteroatom/oxygen not part of this single tellurocyanate
group.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._substituents import name_branch

_YNE_ORDER = 3.0
_TELLURIUM = 52


def _tellurocyanate_cores(mol):
    """List of (tellurium_idx, nitrile_carbon_idx, alkyl_carbon_idx) for
    every R-Te-C#N pattern: a tellurium bonded to exactly one carbon (R) by
    a single bond and one carbon (the nitrile carbon) by a single bond,
    that nitrile carbon triple-bonded to a terminal nitrogen and bonded to
    nothing else."""
    cores = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _TELLURIUM or atom.GetDegree() != 2:
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


def has_tellurocyanate_shape(mol) -> bool:
    return bool(_tellurocyanate_cores(mol))


def name_tellurocyanate(mol) -> str:
    cores = _tellurocyanate_cores(mol)
    if len(cores) != 1:
        raise UnsupportedStructure(
            "exactly one tellurocyanate (-Te-C#N) group is required; zero "
            "or multiple such groups are not supported yet"
        )
    tellurium_idx, nitrile_c_idx, alkyl_c_idx = cores[0]

    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-attached tellurocyanate is out of scope for this "
            "acyclic-only module"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (nitrogen_idx,) = (
        n.GetIdx() for n in mol.GetAtomWithIdx(nitrile_c_idx).GetNeighbors() if n.GetIdx() != tellurium_idx
    )
    excluded = {tellurium_idx, nitrile_c_idx, nitrogen_idx}
    for atom in mol.GetAtoms():
        if atom.GetIdx() in excluded:
            continue
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms other than this single tellurocyanate group "
                "are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0 or atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic rings are out of scope for this module")

    non_single = [b for b in non_single_bonds(mol) if nitrile_c_idx not in (b[0], b[1])]
    if non_single:
        raise UnsupportedStructure("unsaturation in the R group is not supported yet")

    r_name, _ = name_branch(adjacency(mol), alkyl_c_idx, tellurium_idx, {})
    return f"{r_name} tellurocyanate"
