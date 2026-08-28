"""Naming of sulfoxides (R-S(=O)-R'), restricted to two acyclic unbranched
saturated hydrocarbon chains hung off a single sulfinyl sulfur, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-63.6 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  sulfoxide's preferred IUPAC name is formed substitutively (method (1)) by
  prefixing the acyl group name R'-S(=O)- to the parent hydride name for R
  (P-65.3.2.2.2) -- not, as PubChem's auto-generated name for these
  structures might suggest, by treating R' as a plain alkyl substituent with
  an 'oxy'-like 'sulfinyl' suffix tacked directly on. The acyl group name is
  built from the corresponding sulfinic acid name (`_sulfinic_acid.py`'s own
  'methanesulfinic acid' pattern), e.g. R'=ethyl -> 'ethanesulfinyl' (from
  'ethanesulfinic acid'), *not* 'ethylsulfinyl'.
- Confirmed via two of the Blue Book's own worked examples for this exact
  acyclic R-S(=O)-R' shape: '1-(ethanesulfinyl)butane (PIN)' for
  ethyl butyl sulfoxide (R=butyl, R'=ethyl), and (for the analogous sulfone,
  see `_sulfone.py`) '(ethanesulfonyl)ethane (PIN)' for diethyl sulfone.
- The acyl-derived prefix is always cited in enclosing marks, e.g.
  '(ethanesulfinyl)', regardless of chain length or whether a locant
  precedes it -- confirmed by every worked example in P-63.6, including the
  two-carbon symmetric case '(ethanesulfonyl)ethane (PIN)', which has *no*
  locant despite the enclosed prefix (P-14.3.4.2(b): a homogeneous
  two-carbon parent chain bearing a single substituent has only one
  possible structure regardless of numbering direction or how that
  substituent's own name is formed).
- Choice of parent side (R vs R') follows P-44.3, mirroring `_ether.py`/
  `_peroxide.py`: the side with the greater number of skeletal (carbon)
  atoms is the parent hydride; the other becomes the acyl-derived prefix.

This module writes its own minimal locant logic instead of reusing
`_acyclic.py`'s shared `name_from_carbon_graph` machinery, since that
machinery's two-carbon locant-omission shortcut is currently gated on the
substituent *not* being a parenthesized ("compound") one -- a restriction
the '(ethanesulfonyl)ethane (PIN)' example above shows doesn't hold in
general. Fixing that shared helper is out of scope here; this module simply
implements the (very small) locant rule directly for its own single-acyl-
substituent case.

Explicitly out of scope (raise `UnsupportedStructure`):
- A branched R or R' (both sides are required to be a straight,
  unbranched alkyl chain -- unlike `_ether.py`, this isn't just an
  enclosure-mark limitation on one side, since a branched R' would also
  need a locant *inside* its own acid-derived name, e.g.
  'propane-2-sulfinyl', which this module doesn't construct).
- More than one sulfoxide group, any other heteroatom, any unsaturation,
  any ring, or any halogen substituent.
- Selenium/tellurium chalcogen analogues (selenoxide/telluroxide).
"""

from ._common import UnsupportedStructure, non_single_bonds
from ._numerals import alkane_name


def _sulfinyl_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfoxide group: bonded to exactly two
    carbons and one double-bonded (terminal) oxygen (degree 3)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 3:
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


def has_sulfoxide_shape(mol) -> bool:
    return bool(_sulfinyl_sulfur_atoms(mol))


def _unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None if
    the chain branches, rings, or leaves carbon at any point."""
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


def name_sulfoxide(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 16):
            raise UnsupportedStructure(
                "heteroatoms other than the sulfoxide's own sulfur and "
                "oxygen are not supported yet (P-63.6 is restricted to a "
                "plain acyclic sulfoxide here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    sulfurs = _sulfinyl_sulfur_atoms(mol)
    if len(sulfurs) != 1:
        raise UnsupportedStructure("more than one sulfoxide group is out of scope for this module")
    (sulfur,) = sulfurs
    sulfinyl_atom_idxs = {sulfur.GetIdx()} | {n.GetIdx() for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 8}
    if any(
        a not in sulfinyl_atom_idxs and b not in sulfinyl_atom_idxs for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-63.6's scope "
            "here is limited to two saturated chains)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    s_idx = sulfur.GetIdx()
    c1_idx, c2_idx = (n.GetIdx() for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 6)

    len1 = _unbranched_chain_length(mol, c1_idx, s_idx)
    len2 = _unbranched_chain_length(mol, c2_idx, s_idx)
    if len1 is None or len2 is None:
        raise UnsupportedStructure(
            "a branched R or R' group is out of scope for this module (see "
            "module docstring)"
        )

    parent_len, acyl_len = (len1, len2) if len1 >= len2 else (len2, len1)

    acyl_prefix = f"({alkane_name(acyl_len)}sulfinyl)"
    parent = alkane_name(parent_len)
    if parent_len <= 2:
        # P-14.3.4.2(a)/(b): a mononuclear parent's locant is always '1' and
        # never cited; a homogeneous two-carbon parent bearing this single
        # substituent has only one possible structure either way -- see
        # module docstring for why this doesn't depend on the substituent
        # being enclosed.
        return acyl_prefix + parent
    return f"1-{acyl_prefix}{parent}"
