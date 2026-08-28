"""Naming of simple ammonium cations (the '-aminium'/'azanium' suffix
family, P-73.1.1.2), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-73.1.1.2 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf):
  a cation formed by adding a hydron to a parent hydride is named by
  changing the parent hydride name's terminal 'e' to the suffix 'ium'.
  Applied to a primary amine (`_amine.py`'s '-amine' names, P-33.1), this
  gives e.g. 'methanamine' -> 'methanaminium', 'propan-1-amine' ->
  'propan-1-aminium'. Unsubstituted NH4+ is the hydron-added form of
  'azane' (NH3) itself, giving 'azanium' -- both the Blue Book PIN and
  PubChem's own auto-generated name agree on this one (an unusually
  complete match for this project).
- This module reuses `_amine.py`'s existing chain/ring naming wholesale:
  a primary R-NH3+ ammonium is neutralized to the equivalent R-NH2 amine
  (formal charge 0, one fewer hydrogen), named via `name_amine`, and the
  resulting name's terminal 'e' is replaced with 'ium'. No new locant or
  substituent-ordering logic is needed -- `_amine.py`'s existing scope
  (acyclic/monocyclic, halogen coexistence) carries over unchanged.

Explicitly out of scope (raise `UnsupportedStructure`):
- Secondary/tertiary/quaternary ammonium (nitrogen bonded to more than one
  carbon) -- `_amine.py` itself does not support secondary/tertiary amines
  yet, so there is no '-amine' name to derive 'ium' from.
- Any ammonium nitrogen not shaped like a simple, singly-charged R-NH3+ or
  NH4+ (e.g. formal charge other than +1, more than one charged atom,
  isotopic modification, a nitrogen double/triple-bonded to carbon).
- Any other cation-forming parent (oxonium R3O+, sulfonium R3S+, ...) --
  each is a separate P-73 subsection with its own derivation rule.
"""

from rdkit import Chem

from ._amine import name_amine
from ._common import UnsupportedStructure


def has_ammonium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged nitrogen shaped
    like a genuine ammonium (NH4+, or a nitrogen singly bonded to exactly
    one carbon with three hydrogens). Used by `core.py` to route here
    before `_amine.py`, which rejects any charged atom outright.

    Deliberately narrower than "any charged nitrogen exists": a nitro
    group's canonical Lewis structure (-[N+](=O)[O-]) and an isocyanide's
    (-[N+]#[C-]) both also carry a formally charged nitrogen, but neither
    is degree-1/singly-bonded-to-carbon/three-H shaped, so this predicate
    correctly leaves them to their own modules."""
    charged_nitrogens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1
    ]
    if len(charged_nitrogens) != 1:
        return False
    nitrogen = charged_nitrogens[0]
    if nitrogen.GetIsotope() != 0:
        return False
    degree = nitrogen.GetDegree()
    if degree == 0:
        return nitrogen.GetTotalNumHs() == 4
    if degree == 1:
        (neighbor,) = nitrogen.GetNeighbors()
        bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), neighbor.GetIdx())
        return (
            neighbor.GetAtomicNum() == 6
            and bond.GetBondTypeAsDouble() == 1.0
            and nitrogen.GetTotalNumHs() == 3
        )
    return False


def name_ammonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_nitrogens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() != 0
    ]
    (nitrogen,) = charged_nitrogens
    if nitrogen.GetFormalCharge() != 1 or nitrogen.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "ammonium nitrogen is supported (P-73.1.1.2)"
        )

    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() != nitrogen.GetIdx()]
    if not other_atoms:
        if nitrogen.GetDegree() != 0 or nitrogen.GetTotalNumHs() != 4:
            raise UnsupportedStructure("an unsupported unsubstituted ammonium shape")
        return "azanium"

    if nitrogen.GetDegree() != 1 or nitrogen.GetTotalNumHs() != 3:
        raise UnsupportedStructure(
            "only an isolated primary ammonium (R-NH3+, nitrogen bonded to "
            "exactly one carbon, three hydrogens) is supported yet; "
            "secondary/tertiary/quaternary ammonium is out of scope "
            "(P-73.1.1.2)"
        )
    (neighbor,) = nitrogen.GetNeighbors()
    if neighbor.GetAtomicNum() != 6:
        raise UnsupportedStructure("an ammonium nitrogen must be attached to a carbon atom")
    bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), neighbor.GetIdx())
    if bond.GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("the ammonium nitrogen must be singly bonded to carbon")

    neutral_rw = Chem.RWMol(mol)
    neutral_nitrogen = neutral_rw.GetAtomWithIdx(nitrogen.GetIdx())
    neutral_nitrogen.SetFormalCharge(0)
    neutral_nitrogen.SetNoImplicit(True)
    neutral_nitrogen.SetNumExplicitHs(2)
    neutral_mol = neutral_rw.GetMol()
    Chem.SanitizeMol(neutral_mol)

    amine_name = name_amine(neutral_mol)
    return amine_name[:-1] + "ium"
