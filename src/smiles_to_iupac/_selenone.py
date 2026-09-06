"""Naming of selenones (R-Se(=O)(=O)-R'), the selenium analogue of a
sulfone, restricted to two acyclic unbranched saturated hydrocarbon chains
hung off a single selenonyl selenium, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-63.6 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  selenone's preferred IUPAC name is formed the same way as a selenoxide's
  (see `_selenoxide.py`'s module docstring, which this mirrors exactly) --
  substitutively, prefixing the acyl group R'-Se(=O)(=O)- (built from the
  corresponding selenonic acid name, e.g. R'=ethyl -> 'ethaneselenonyl') to
  the parent hydride name for R.
- Confirmed via the Blue Book's own worked example for this exact acyclic
  R-Se(=O)(=O)-R' shape: '7-(benzeneselenonyl)quinoline (PIN)'.

Explicitly out of scope (raise `UnsupportedStructure`): same list as
`_selenoxide.py` -- a branched R or R', more than one selenone group, any
other heteroatom, any unsaturation, any ring, any halogen substituent, and
the tellurium analogue (tellurone).
"""

from rdkit import Chem

from ._common import UnsupportedStructure, non_single_bonds
from ._numerals import alkane_name


def _unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None if
    the chain branches, rings, or leaves carbon at any point. Mirrors
    `_selenoxide.py`'s identical helper."""
    length = 0
    previous = exclude_idx
    current = root_idx
    while True:
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return None
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous]
        length += 1
        if not neighbors:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]


def _selenonyl_selenium_atoms(mol):
    """Selenium atoms shaped like a selenone group: bonded to exactly two
    carbons and two double-bonded (terminal) oxygens (degree 4)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 34 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 2 or len(oxygens) != 2:
            continue
        if any(
            mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() != 2.0 or o.GetDegree() != 1
            for o in oxygens
        ):
            continue
        matches.append(atom)
    return matches


def has_selenone_shape(mol) -> bool:
    return bool(_selenonyl_selenium_atoms(mol))


def name_selenone(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 34):
            raise UnsupportedStructure(
                "heteroatoms other than the selenone's own selenium and "
                "oxygens are not supported yet (P-63.6 is restricted to a "
                "plain acyclic selenone here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    seleniums = _selenonyl_selenium_atoms(mol)
    if len(seleniums) != 1:
        raise UnsupportedStructure("more than one selenone group is out of scope for this module")
    (selenium,) = seleniums
    selenonyl_atom_idxs = {selenium.GetIdx()} | {n.GetIdx() for n in selenium.GetNeighbors() if n.GetAtomicNum() == 8}
    if any(
        a not in selenonyl_atom_idxs and b not in selenonyl_atom_idxs for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-63.6's scope "
            "here is limited to two saturated chains)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    se_idx = selenium.GetIdx()
    c1_idx, c2_idx = (n.GetIdx() for n in selenium.GetNeighbors() if n.GetAtomicNum() == 6)

    len1 = _unbranched_chain_length(mol, c1_idx, se_idx)
    len2 = _unbranched_chain_length(mol, c2_idx, se_idx)
    if len1 is None or len2 is None:
        raise UnsupportedStructure(
            "a branched R or R' group is out of scope for this module (see "
            "module docstring)"
        )

    parent_len, acyl_len = (len1, len2) if len1 >= len2 else (len2, len1)

    acyl_prefix = f"({alkane_name(acyl_len)}selenonyl)"
    parent = alkane_name(parent_len)
    if parent_len <= 2:
        return acyl_prefix + parent
    return f"1-{acyl_prefix}{parent}"
