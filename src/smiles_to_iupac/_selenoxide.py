"""Naming of selenoxides (R-Se(=O)-R'), the selenium analogue of a sulfoxide,
restricted to two acyclic unbranched saturated hydrocarbon chains hung off a
single seleninyl selenium, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-63.6 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  "Selenium and tellurium analogues are named in the same way using acyl
  groups derived from the appropriate seleninic, selenonic, tellurinic, and
  telluronic acids" -- i.e. exactly `_sulfoxide.py`'s substitutive rule
  (prefixing the acyl group R'-Se(=O)- to the parent hydride name for R),
  with 'seleninyl' (from 'seleninic acid') standing in for 'sulfinyl'.
- Confirmed via the Blue Book's own worked example for this exact acyclic
  R-Se(=O)-R' shape: '(ethaneseleninyl)benzene (PIN)' for ethyl phenyl
  selenoxide.
- Choice of parent side and the two-carbon no-locant rule are identical to
  `_sulfoxide.py` -- see that module's docstring for the P-44.3/P-14.3.4.2(b)
  reasoning, which applies unchanged here.

This module mirrors `_sulfoxide.py` structurally with Se in place of S;
see that module for the rationale behind writing its own minimal locant
logic instead of reusing `_acyclic.py`.

Explicitly out of scope (raise `UnsupportedStructure`): same list as
`_sulfoxide.py` -- a branched R or R', more than one selenoxide group, any
other heteroatom, any unsaturation, any ring, and any halogen substituent.
(The tellurium analogue, telluroxide, is handled separately by
`_telluroxide.py`.)
"""

from rdkit import Chem

from ._common import UnsupportedStructure, non_single_bonds, unbranched_chain_length
from ._numerals import alkane_name


def _seleninyl_selenium_atoms(mol):
    """Selenium atoms shaped like a selenoxide group: bonded to exactly two
    carbons and one double-bonded (terminal) oxygen (degree 3)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 34 or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 2 or len(oxygens) != 1:
            continue
        (oxygen,) = oxygens
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx())
        if bond.GetBondTypeAsDouble() != 2.0 or oxygen.GetDegree() != 1:
            continue
        matches.append(atom)
    return matches


def has_selenoxide_shape(mol) -> bool:
    return bool(_seleninyl_selenium_atoms(mol))


def name_selenoxide(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 34):
            raise UnsupportedStructure(
                "heteroatoms other than the selenoxide's own selenium and "
                "oxygen are not supported yet (P-63.6 is restricted to a "
                "plain acyclic selenoxide here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    seleniums = _seleninyl_selenium_atoms(mol)
    if len(seleniums) != 1:
        raise UnsupportedStructure("more than one selenoxide group is out of scope for this module")
    (selenium,) = seleniums
    seleninyl_atom_idxs = {selenium.GetIdx()} | {n.GetIdx() for n in selenium.GetNeighbors() if n.GetAtomicNum() == 8}
    if any(
        a not in seleninyl_atom_idxs and b not in seleninyl_atom_idxs for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-63.6's scope "
            "here is limited to two saturated chains)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    se_idx = selenium.GetIdx()
    c1_idx, c2_idx = (n.GetIdx() for n in selenium.GetNeighbors() if n.GetAtomicNum() == 6)

    len1 = unbranched_chain_length(mol, c1_idx, se_idx)
    len2 = unbranched_chain_length(mol, c2_idx, se_idx)
    if len1 is None or len2 is None:
        raise UnsupportedStructure(
            "a branched R or R' group is out of scope for this module (see "
            "module docstring)"
        )

    parent_len, acyl_len = (len1, len2) if len1 >= len2 else (len2, len1)

    acyl_prefix = f"({alkane_name(acyl_len)}seleninyl)"
    parent = alkane_name(parent_len)
    if parent_len <= 2:
        return acyl_prefix + parent
    return f"1-{acyl_prefix}{parent}"
