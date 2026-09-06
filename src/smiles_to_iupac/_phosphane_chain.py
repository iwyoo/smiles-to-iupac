"""Naming of unbranched, unsubstituted acyclic phosphorus hydride chains
(phosphane, diphosphane, triphosphane, ...), per the IUPAC 2013
Recommendations ("the Blue Book"):

- Chapter P-6a (https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf, line 7709):
  a chain of n phosphorus atoms (all single P-P bonds, each terminal P
  bearing two H, each internal P bearing one H -- i.e. PnH(n+2), the
  phosphorus analogue of `_silane_chain.py`'s SinH(2n+2) silane chain) is
  named with the P-14.2.1 numerical-term multiplying prefix ('di', 'tri',
  ...) plus 'phosphane'; n=1 (PH3) is simply 'phosphane', already handled
  by `_phosphane.py`'s own mononuclear case. Confirmed directly from the
  primary text's own worked example: "diphosphane (preselected name)".
- **PubChem gives a non-PIN name for this exact case**: `PP` -> PubChem
  CID 139283 "phosphanylphosphane", `PPP` -> CID 139510
  "bis(phosphanyl)phosphane" -- neither matches the Blue Book's own
  multiplying-prefix pattern. This module follows the primary text
  instead, the same established practice `_silane_chain.py` already set
  for its own "disilane"/"trisilane" pattern (this project has repeatedly
  found PubChem's auto-generated names disagreeing with the actual PIN
  elsewhere too).
- Boron (B-B) chains are deliberately out of scope here: diborane uses a
  structurally different naming system ('diborane(4)'/'diborane(6)',
  P6a.pdf lines 5472-5491) driven by non-classical bridging-hydrogen
  bonding, unlike phosphorus's plain classical-valence chain -- a
  separate, unverified problem.

Explicitly out of scope (raise `UnsupportedStructure`): any branching, any
substituent (halogen included), any atom other than phosphorus, any bond
order other than single, and any ring.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._numerals import numerical_term

_PHOSPHORUS = 15


def has_phosphane_chain_shape(mol) -> bool:
    """True iff every atom in `mol` is phosphorus (the shape this module
    accepts -- everything else is left to `core.py`'s other branches)."""
    atoms = mol.GetAtoms()
    return len(atoms) > 0 and all(atom.GetAtomicNum() == _PHOSPHORUS for atom in atoms)


def name_phosphane_chain(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure(
                "a non-single P-P bond is not supported yet (see P-21.2.1)"
            )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("cyclic phosphorus skeletons are not supported yet (see P-22)")

    graph = adjacency(mol)
    n = len(graph)
    if n == 1:
        return "phosphane"
    degrees = sorted(len(neighbors) for neighbors in graph.values())
    if degrees != [1, 1] + [2] * (n - 2):
        raise UnsupportedStructure(
            "a branched phosphorus chain is not supported yet (see "
            "P-21.2.1's general branched-chain substitutive naming, not "
            "yet implemented for phosphanes)"
        )
    # 'phosphane' starts with a consonant, so the numerical term's
    # terminal 'a' is kept unchanged, same as `_silane_chain.py`'s
    # 'tetra' + 'silane' -> 'tetrasilane'.
    return numerical_term(n) + "phosphane"
