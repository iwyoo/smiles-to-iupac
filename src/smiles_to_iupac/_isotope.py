"""Naming of isotopically substituted methane (the '(<nuclide>)' isotope
descriptor, P-82.2.1), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-82.2.1 (Chapter P-8, https://iupac.qmul.ac.uk/BlueBook/P8.html): the name
  of an isotopically substituted compound is formed by adding the nuclide
  symbol(s), enclosed in parentheses, preceded by any necessary locant(s),
  before the part of the name that is isotopically substituted. Confirmed
  worked examples: CH3-2H -> '(2H1)methane (PIN)'; (14C)methane (PIN);
  trichloro(12C)methane (PIN); dichloro(2H2)methane (PIN). The isotope's
  count is always cited as a right subscript to its symbol, even for a
  single atom ('2H1', never a bare '2H') -- "the number of atoms
  substituted is always specified as a right subscript ... even in case of
  monosubstitution".
- This module covers only a mononuclear methane parent: no locant is ever
  needed here, the same way `_acyclic.py`'s P-14.3.4.2(a) never cites a
  locant on a single-carbon parent (e.g. 'chloromethane'). Longer chains
  (ethane, propane, ...) would additionally require confirming exactly when
  a locant is cited for a symmetric *unmodified* parent name -- P-82.2's own
  worked example for a multi-carbon chain, '(2-13C)ethan-1-ol', only covers
  an already-locanted parent ('ethan-1-ol', not plain 'ethane'), so it
  doesn't settle the unlabeled-parent case -- out of scope until a second
  worked example confirms it.
- Within methane, two independent isotope axes are supported: 1-4 deuterium
  (2H) atoms replacing hydrogen, or the skeletal carbon itself being 12C/
  13C/14C -- each optionally coexisting with the already-supported halogen
  substituent prefixes (P-35.2.1).

Explicitly out of scope (raise `UnsupportedStructure`):
- Any parent other than methane (CH4) itself -- see above.
- Any isotope other than deuterium (2H) among hydrogen atoms (e.g. tritium).
- The carbon isotope and deuterium substitution occurring at the same time
  (P-82.3's rules for combining more than one nuclide symbol in one name),
  or any isotopically labeled halogen.
- Any heteroatom (other than halogens), charge, or additional isotopic
  modification.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._numerals import numerical_term

_CARBON_ISOTOPES = {12, 13, 14}
_MAX_DEUTERIUMS = 4


def has_isotope_shape(mol) -> bool:
    """True if the molecule contains any isotopically labeled atom at all
    (deuterium, a heavy carbon isotope, or otherwise), regardless of
    whether the rest of the molecule is in scope. Used by `core.py` to
    route here before every other branch, none of which recognize an
    explicit, isotopically labeled atom (RDKit represents 2H as its own
    atom, atomic number 1) at all."""
    return any(atom.GetIsotope() != 0 for atom in mol.GetAtoms())


def name_isotope(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) != 1:
        raise UnsupportedStructure(
            "this module only supports a single mononuclear methane parent; "
            "other parent hydrides are not supported yet (P-82.2.1)"
        )
    carbon = carbons[0]
    if carbon.GetFormalCharge() != 0 or carbon.GetIsAromatic():
        raise UnsupportedStructure("an unsupported carbon shape for methane")

    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() != carbon.GetIdx()]
    if any(atom.GetAtomicNum() not in {1, *HALOGEN_PREFIXES} for atom in other_atoms):
        raise UnsupportedStructure(
            "this module only supports deuterium and halogen substituents "
            "on methane; other atoms are not supported yet (P-82.2.1)"
        )
    for atom in other_atoms:
        bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
        if bond is None or bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("substituents on methane must be singly bonded")
        if atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetIsotope() != 0:
            raise UnsupportedStructure("an isotopically labeled halogen is not supported yet")

    hydrogens = [atom for atom in other_atoms if atom.GetAtomicNum() == 1]
    non_deuteriums = [atom for atom in hydrogens if atom.GetIsotope() != 2]
    if non_deuteriums:
        raise UnsupportedStructure(
            "only deuterium (2H) is supported among hydrogen isotopes; "
            "other hydrogen isotopes (e.g. tritium) are not supported yet"
        )
    deuterium_count = len(hydrogens)

    carbon_isotope = carbon.GetIsotope()
    if carbon_isotope != 0 and carbon_isotope not in _CARBON_ISOTOPES:
        raise UnsupportedStructure(
            "only 12C/13C/14C skeletal carbon isotopes are supported (P-82.2.1)"
        )
    if carbon_isotope != 0 and deuterium_count:
        raise UnsupportedStructure(
            "a skeletal carbon isotope combined with deuterium substitution "
            "is not supported yet (P-82.3)"
        )
    if deuterium_count > _MAX_DEUTERIUMS:
        raise UnsupportedStructure("methane cannot carry more than 4 deuterium atoms")

    if carbon_isotope != 0:
        isotope_descriptor = f"{carbon_isotope}C"
    elif deuterium_count:
        isotope_descriptor = f"2H{deuterium_count}"
    else:
        raise UnsupportedStructure("no isotopically labeled atom found")

    halogens = [atom for atom in other_atoms if atom.GetAtomicNum() in HALOGEN_PREFIXES]
    halogen_counts: dict[str, int] = {}
    for atom in halogens:
        prefix = HALOGEN_PREFIXES[atom.GetAtomicNum()]
        halogen_counts[prefix] = halogen_counts.get(prefix, 0) + 1
    halogen_prefix = "".join(
        prefix if count == 1 else numerical_term(count) + prefix
        for prefix, count in sorted(halogen_counts.items())
    )

    return f"{halogen_prefix}({isotope_descriptor})methane"
