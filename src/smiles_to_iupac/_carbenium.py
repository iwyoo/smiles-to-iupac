"""Naming of simple carbenium cations ('methylium', 'propylium',
'cyclobutylium', ...), per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-73.2.2.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf),
  the "specific method": "A cation formally derived by the removal of a
  hydride ion, H-, ... from a terminal atom of a saturated unbranched
  acyclic hydrocarbon [or] a saturated monocyclic hydrocarbon ... is named
  by replacing the 'ane' ending in the name of the parent hydride by the
  suffix 'ylium'." Confirmed worked examples (`tmp/bluebook/P7.txt`):
  `[CH3+]` -> 'methylium (PIN)'; a terminal cation on propane ->
  'propylium (PIN)'; a cyclobutane ring cation -> 'cyclobutylium (PIN)'.
  This is structurally identical to `_radical.py`'s own "ane"->"yl"
  mechanism (P-71.2.1.1), just with 'ylium' instead of 'yl' -- this
  module reuses `_numerals.py`'s `alkyl_name` the same way, appending
  'ium' to its own output.

Explicitly out of scope (raise `UnsupportedStructure`), mirroring
`_radical.py`'s own original scope before its later branch-point
extension:
- A cation carbon that is itself a branch point (not a chain terminus) --
  P-73.2.2.1.2's "general method" names this by citing the cation's own
  locant directly on the full (possibly branched) parent hydride name
  (e.g. 'heptamethyltrisilan-2-ylium (PIN)'), which would need
  `_substituents.py`'s `name_branch` for any branched substituent --
  currently pre-PIN (P-29 blocker, see the roadmap) -- so this shape is
  deferred to a follow-up task, exactly as `_radical.py`'s own P-29.3.2.2
  branch-point handling was.
- More than one cationic center, a cation on a polycyclic/spiro skeleton
  or an aromatic ring, or coexisting with any heteroatom, halogen,
  unsaturation, or isotopic modification.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._numerals import alkyl_name


def has_carbenium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged carbon shaped
    like a genuine carbenium (a carbon with total substituent+hydrogen
    count of 3, i.e. one fewer bond than neutral carbon's four). Used by
    `core.py` to route here before any other branch, none of which
    recognize a charged atom."""
    charged_carbons = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 1
    ]
    if len(charged_carbons) != 1:
        return False
    carbon = charged_carbons[0]
    if carbon.GetIsotope() != 0:
        return False
    degree = carbon.GetDegree()
    if degree > 3 or carbon.GetTotalNumHs() + degree != 3:
        return False
    return all(
        n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in carbon.GetNeighbors()
    )


def name_carbenium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_carbons = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() != 0
    ]
    (cation,) = charged_carbons
    if cation.GetFormalCharge() != 1 or cation.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "carbenium carbon is supported (P-73.2.2.1.1)"
        )

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "only an all-carbon skeleton is supported yet (P-73.2.2.1.1)"
            )
        if atom.GetIdx() != cation.GetIdx() and atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("more than one charged atom is not supported yet")
        if atom.GetIsotope() != 0:
            raise UnsupportedStructure("isotopically modified atoms are not supported yet")
        if atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic rings are out of scope for this module")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("unsaturated skeletons are not supported yet")

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_chain_carbenium(mol, cation)
    if num_rings == 1:
        return _name_ring_carbenium(mol, ring_info)
    raise UnsupportedStructure("polycyclic and spiro carbenium cations are not supported yet")


def _name_chain_carbenium(mol, cation) -> str:
    if mol.GetNumAtoms() == 1:
        # P-73.2.2.1.1's own wording covers this directly: a mononuclear
        # parent hydride (methane) -- methylium (CH3+) has no bond to
        # another atom, so degree 0 rather than 1.
        return alkyl_name(1) + "ium"
    if cation.GetDegree() != 1:
        raise UnsupportedStructure(
            "a carbenium carbon that is itself a branch point is out of "
            "scope for this module (P-73.2.2.1.2, the 'general method')"
        )
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            raise UnsupportedStructure(
                "a branched chain is out of scope for this module "
                "(P-73.2.2.1.2, the 'general method')"
            )
    return alkyl_name(mol.GetNumAtoms()) + "ium"


def _name_ring_carbenium(mol, ring_info) -> str:
    (ring_atoms,) = ring_info.AtomRings()
    if len(ring_atoms) != mol.GetNumAtoms():
        raise UnsupportedStructure(
            "a substituent hanging off the ring (other than the "
            "carbenium cation itself) is out of scope for this module"
        )
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 2:
            raise UnsupportedStructure(
                "a monocyclic carbenium ring must otherwise be "
                "unsubstituted (P-73.2.2.1.1)"
            )
    return "cyclo" + alkyl_name(len(ring_atoms)) + "ium"
