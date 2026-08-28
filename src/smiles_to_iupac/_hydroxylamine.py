"""Naming of hydroxylamine (H2N-OH) and its O-substituted derivatives, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-68.3.1.1.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'hydroxylamine' is a preselected (retained) name for H2N-OH, and is a
  functional parent compound that allows substitution on its oxygen atom.
- P-68.3.1.1.1.2: substitution on the oxygen atom only (H2N-O-R) is
  expressed as 'O-<R>hydroxylamine' -- confirmed via a worked example,
  'O-methylhydroxylamine (PIN)' for H2N-O-CH3.

Substitution on the *nitrogen* atom (R-NH-OH, RR'N-OH) is out of scope
here: P-68.3.1.1.1.1 gives its PIN as an N-hydroxy derivative of the
corresponding amine (e.g. CH3-NH-OH is 'N-hydroxymethanamine (PIN)', not
'N-methylhydroxylamine') -- a different parent-selection mechanism
(composing with `_amine.py`'s chain-naming, not this module's mononuclear
hydroxylamine parent) that this module doesn't implement. N,O-disubstitution
(P-68.3.1.1.1.3) is likewise out of scope, since its PIN is also an
N-substituted amine (e.g. CH3-NH-O-CH3 is 'N-methoxymethanamine (PIN)').
Molecules shaped like either therefore don't match `has_hydroxylamine_shape`
and fall through to the existing amine module's own rejection.

Explicitly out of scope (raise `UnsupportedStructure` from the shape check
returning False, or fall through to another module):
- Any N-substitution (nitrogen degree != 1) -- different PIN mechanism, see
  above.
- A branched (compound) O-substituent -- mirrors `_ether.py`'s scope limit.
- Any unsaturation, any ring, or any heteroatom other than the hydroxylamine's
  own N and O.
"""

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._substituents import name_branch


def _hydroxylamine_atoms(mol):
    """(n, o) atoms of a plain hydroxylamine skeleton -- one nitrogen, one
    oxygen, singly bonded to each other, nitrogen unsubstituted (degree 1,
    bonded only to the oxygen) -- or None if `mol` isn't shaped this way."""
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(nitrogens) != 1 or len(oxygens) != 1:
        return None
    n, o = nitrogens[0], oxygens[0]
    bond = mol.GetBondBetweenAtoms(n.GetIdx(), o.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return None
    if n.GetDegree() != 1:
        return None
    if o.GetDegree() not in (1, 2):
        return None
    return n, o


def has_hydroxylamine_shape(mol) -> bool:
    return _hydroxylamine_atoms(mol) is not None


def name_hydroxylamine(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 7, 8):
            raise UnsupportedStructure(
                "heteroatoms other than hydroxylamine's own nitrogen and "
                "oxygen are not supported yet (P-68.3.1.1.1 is restricted "
                "to a plain hydroxylamine skeleton here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-68.3.1.1.1's "
            "O-substituent scope here is limited to a saturated chain)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    n, o = _hydroxylamine_atoms(mol)
    if o.GetDegree() == 1:
        return "hydroxylamine"

    n_idx, o_idx = n.GetIdx(), o.GetIdx()
    r_root = next(nbr.GetIdx() for nbr in o.GetNeighbors() if nbr.GetIdx() != n_idx)
    graph = adjacency(mol)
    r_name, r_compound = name_branch(graph, r_root, o_idx, {})
    if r_compound:
        raise UnsupportedStructure(
            "a branched O-substituent is not supported yet (mirrors "
            "_ether.py's own enclosure-mark limitation)"
        )
    return f"O-{r_name}hydroxylamine"
