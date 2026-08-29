"""Naming of simple phosphonium cations (the '-phosphanium' suffix,
P-73.1.1.2), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-73.1.1.2 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf):
  a cation formed by adding a hydron to a parent hydride is named by
  changing the parent hydride name's terminal 'e' to the suffix 'ium' --
  the same rule `_ammonium.py` already applies to amines, just for
  phosphorus's own parent hydride 'phosphane' (PH3) instead of nitrogen's
  'azane'. Since `_phosphane.py`'s own substitutive naming already covers
  0-3 plain alkyl substituents on phosphorus (unlike `_amine.py`, which
  only covers a single primary substituent), this module reuses that
  scope wholesale: a phosphonium cation is neutralized to the equivalent
  phosphane (formal charge 0, one fewer hydrogen), named via
  `name_simple_phosphane`, and the resulting name's terminal 'e' is
  replaced with 'ium'.
- Confirmed via PubChem structure match: `[PH4+]` -> "phosphanium",
  `C[PH3+]` -> "methylphosphanium", `C[PH2+]C` -> "dimethylphosphanium",
  `C[PH+](C)C` -> "trimethylphosphanium". A fourth substituent
  (`C[P+](C)(C)C` -> "tetramethylphosphanium", a genuine quaternary
  phosphonium salt) is structurally confirmed too but is NOT reachable by
  this module's neutralize-then-rename approach -- a neutral phosphorus
  atom cannot carry four substituents at all, so there is no phosphane
  name to derive it from; naming that shape needs its own substitutive
  logic and is deferred to a follow-up task.

Explicitly out of scope (raise `UnsupportedStructure`):
- Quaternary phosphonium (phosphorus bonded to four carbons) -- see above.
- Any phosphonium phosphorus not shaped like PH4+ or a phosphorus bonded
  to 1-3 carbons (with the remaining valence as hydrogens) -- e.g. formal
  charge other than +1, more than one charged atom, isotopic
  modification, a halogen or other heteroatom substituent, a
  branched/unsaturated/aromatic/ring-bearing substituent (inherited
  unchanged from `_phosphane.py`'s own scope, since this module's
  validation is entirely delegated to it after neutralization).
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._phosphane import name_simple_phosphane


def has_phosphonium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged phosphorus
    shaped like a genuine phosphonium (PH4+, or a phosphorus singly bonded
    to 1-3 carbons with the rest hydrogens). Used by `core.py` to route
    here before `_phosphane.py`, which rejects any charged atom outright."""
    charged_phosphorus = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() == 1
    ]
    if len(charged_phosphorus) != 1:
        return False
    phosphorus = charged_phosphorus[0]
    if phosphorus.GetIsotope() != 0:
        return False
    degree = phosphorus.GetDegree()
    if degree > 3 or phosphorus.GetTotalNumHs() + degree != 4:
        return False
    return all(
        n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in phosphorus.GetNeighbors()
    )


def name_phosphonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_phosphorus = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() != 0
    ]
    (phosphorus,) = charged_phosphorus
    if phosphorus.GetFormalCharge() != 1 or phosphorus.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "phosphonium phosphorus is supported (P-73.1.1.2)"
        )
    degree = phosphorus.GetDegree()
    if degree > 3:
        raise UnsupportedStructure(
            "a quaternary phosphonium (phosphorus bonded to four carbons) "
            "is out of scope for this module -- it has no neutral "
            "phosphane counterpart to derive its name from (P-73.1.1.2)"
        )
    if any(n.GetAtomicNum() != 6 for n in phosphorus.GetNeighbors()):
        raise UnsupportedStructure(
            "a phosphonium substituent other than carbon is out of scope "
            "for this module"
        )
    if any(
        mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
        for n in phosphorus.GetNeighbors()
    ):
        raise UnsupportedStructure("the phosphonium phosphorus must be singly bonded to each substituent")

    neutral_rw = Chem.RWMol(mol)
    neutral_phosphorus = neutral_rw.GetAtomWithIdx(phosphorus.GetIdx())
    neutral_phosphorus.SetFormalCharge(0)
    neutral_phosphorus.SetNoImplicit(True)
    neutral_phosphorus.SetNumExplicitHs(3 - degree)
    neutral_mol = neutral_rw.GetMol()
    Chem.SanitizeMol(neutral_mol)

    phosphane_name = name_simple_phosphane(neutral_mol)
    return phosphane_name[:-1] + "ium"
