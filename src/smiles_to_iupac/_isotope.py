"""Naming of isotopically substituted deuterium-methane (the '(2H1)' isotope
descriptor, P-82.2.1), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-82.2.1 (Chapter P-8, https://iupac.qmul.ac.uk/BlueBook/P8.html): the name
  of an isotopically substituted compound is formed by adding the nuclide
  symbol(s), enclosed in parentheses, preceded by any necessary locant(s),
  before the part of the name that is isotopically substituted. Confirmed
  worked example: CH3-2H -> '(2H1)methane (PIN)'. The isotope's count is
  always cited as a right subscript to its symbol, even for a single atom
  ('2H1', never a bare '2H') -- "the number of atoms substituted is always
  specified as a right subscript ... even in case of monosubstitution".
- This first pass covers only a single deuterium (2H) atom on methane (a
  mononuclear parent): no locant is ever needed here, the same way
  `_acyclic.py`'s P-14.3.4.2(a) never cites a locant on a single-carbon
  parent (e.g. 'chloromethane'). Longer chains (ethane, propane, ...) would
  additionally require confirming exactly when a locant is cited for a
  symmetric *unmodified* parent name -- P-82.2's own worked example for a
  multi-carbon chain, '(2-13C)ethan-1-ol', only covers an already-locanted
  parent ('ethan-1-ol', not plain 'ethane'), so it doesn't settle the
  unlabeled-parent case -- out of scope until a second worked example
  confirms it.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one isotopically labeled atom, or any isotope other than
  deuterium (2H) -- e.g. 13C (a structurally different case involving a
  skeletal atom rather than a substituent atom) is not attempted here.
- Any parent other than methane (CH4) itself -- see above.
- Any heteroatom, charge, or additional isotopic modification.
"""

from rdkit import Chem

from ._common import UnsupportedStructure


def has_isotope_shape(mol) -> bool:
    """True if the molecule contains any isotopically labeled atom at all
    (deuterium or otherwise), regardless of whether the rest of the
    molecule is in scope. Used by `core.py` to route here before every
    other branch, none of which recognize an explicit, isotopically
    labeled hydrogen atom (RDKit represents 2H as its own atom, atomic
    number 1) at all."""
    return any(atom.GetIsotope() != 0 for atom in mol.GetAtoms())


def name_isotope(mol) -> str:
    labeled = [atom for atom in mol.GetAtoms() if atom.GetIsotope() != 0]
    if len(labeled) != 1 or labeled[0].GetAtomicNum() != 1 or labeled[0].GetIsotope() != 2:
        raise UnsupportedStructure(
            "only a single deuterium (2H) atom is supported; zero, "
            "multiple, or non-deuterium isotopic labels are not supported "
            "yet (P-82.2.1)"
        )
    deuterium = labeled[0]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() != deuterium.GetIdx()]
    if len(other_atoms) != 1 or other_atoms[0].GetAtomicNum() != 6:
        raise UnsupportedStructure(
            "this module only supports a single deuterium substituent on "
            "methane; other parent hydrides are not supported yet "
            "(P-82.2.1)"
        )
    carbon = other_atoms[0]
    if carbon.GetFormalCharge() != 0 or carbon.GetIsAromatic():
        raise UnsupportedStructure("an unsupported carbon shape for methane-d1")
    bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), deuterium.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("the deuterium atom must be singly bonded to the carbon")

    return "(2H1)methane"
