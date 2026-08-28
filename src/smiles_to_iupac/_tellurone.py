"""Naming of tellurones (R-Te(=O)(=O)-R'), the tellurium analogue of a
sulfone/selenone, restricted to two acyclic unbranched saturated
hydrocarbon chains hung off a single telluronyl tellurium, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-63.6 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  tellurone's preferred IUPAC name is formed the same way as a
  telluroxide's (see `_telluroxide.py`'s module docstring, which this
  mirrors exactly) -- substitutively, prefixing the acyl group
  R'-Te(=O)(=O)- (built from the corresponding telluronic acid name, e.g.
  R'=ethyl -> 'ethanetelluronyl') to the parent hydride name for R.
- 'telluronyl' is confirmed as a preselected prefix directly in the Blue
  Book text (P6a.txt). No worked example exists for this exact
  R-Te(=O)(=O)-R' shape, and the asymmetric-chain case has no
  PubChem-registered structure (CID 0) -- see `_telluroxide.py`'s docstring
  for why this is a reviewed, not directly verified, mechanical extension.
  The symmetric mononuclear case (dimethyl tellurone) IS PubChem-registered
  (CID 59167377) and is covered by this module's own tests.

Explicitly out of scope (raise `UnsupportedStructure`): same list as
`_telluroxide.py` -- a branched R or R', more than one tellurone group, any
other heteroatom, any unsaturation, any ring, any halogen substituent.
"""

from ._common import UnsupportedStructure, non_single_bonds
from ._numerals import alkane_name


def _unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None if
    the chain branches, rings, or leaves carbon at any point. Mirrors
    `_telluroxide.py`'s identical helper."""
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


def _telluronyl_tellurium_atoms(mol):
    """Tellurium atoms shaped like a tellurone group: bonded to exactly two
    carbons and two double-bonded (terminal) oxygens (degree 4)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 52 or atom.GetDegree() != 4:
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


def has_tellurone_shape(mol) -> bool:
    return bool(_telluronyl_tellurium_atoms(mol))


def name_tellurone(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 52):
            raise UnsupportedStructure(
                "heteroatoms other than the tellurone's own tellurium and "
                "oxygens are not supported yet (P-63.6 is restricted to a "
                "plain acyclic tellurone here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    telluriums = _telluronyl_tellurium_atoms(mol)
    if len(telluriums) != 1:
        raise UnsupportedStructure("more than one tellurone group is out of scope for this module")
    (tellurium,) = telluriums
    telluronyl_atom_idxs = {tellurium.GetIdx()} | {
        n.GetIdx() for n in tellurium.GetNeighbors() if n.GetAtomicNum() == 8
    }
    if any(
        a not in telluronyl_atom_idxs and b not in telluronyl_atom_idxs for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-63.6's scope "
            "here is limited to two saturated chains)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    te_idx = tellurium.GetIdx()
    c1_idx, c2_idx = (n.GetIdx() for n in tellurium.GetNeighbors() if n.GetAtomicNum() == 6)

    len1 = _unbranched_chain_length(mol, c1_idx, te_idx)
    len2 = _unbranched_chain_length(mol, c2_idx, te_idx)
    if len1 is None or len2 is None:
        raise UnsupportedStructure(
            "a branched R or R' group is out of scope for this module (see "
            "module docstring)"
        )

    parent_len, acyl_len = (len1, len2) if len1 >= len2 else (len2, len1)

    acyl_prefix = f"({alkane_name(acyl_len)}telluronyl)"
    parent = alkane_name(parent_len)
    if parent_len <= 2:
        return acyl_prefix + parent
    return f"1-{acyl_prefix}{parent}"
