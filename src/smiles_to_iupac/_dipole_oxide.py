"""Naming of nitrone (imine N-oxide) and nitrile-oxide 1,3-dipoles via
functional-class nomenclature, per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-74.2.1.2 method (2) (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/
  PDF/P7.pdf): a nitrone (R2C=N+(R')-O-, R' possibly H) is named as the
  'N-oxide' of the parent imine it would be if the dipole's own O- were
  removed and the nitrogen's charge neutralized -- "Method (2) leads to
  preferred IUPAC names when one amine oxide is present" (~2920-2929),
  the same construction `_amine_oxide.py` already uses for a plain
  R3N-O amine oxide, applied here to an imine nitrogen instead. Confirmed
  worked example `N-methylethanimine oxide`-shaped (PubChem CID
  59937459), rendered here with the primary source's own N-locant
  convention, `N-methylethanimine N-oxide`, matching the closely related
  `N,N-dimethylmethanamine N-oxide (PIN)` amine-oxide format.
- P-74.2.2.2.1.2 method (2) (~3183-3192): a nitrile oxide (R-C#N+-O-) is
  named the same way -- the parent nitrile's name (with the dipole's O-
  removed and the nitrogen's charge neutralized) plus ' oxide', no
  N-locant needed since a nitrile has only one nitrogen at all. Confirmed
  PIN worked example `acetonitrile oxide (PIN)`, the Blue Book and the Blue Book -- per P-66.5.1.2.1,
  'acetonitrile'/'formonitrile' are themselves the retained PINs for the
  n=2/n=1 parent nitrile ('acetonitrile' comes from the retained-stem table
  in `_common.py`, 'formonitrile' is overridden locally).

Both dipole subtypes share one mechanism: strip the dipole's own oxygen
atom from the molecule, reset the nitrogen's formal charge to neutral (so
the result is exactly the molecule `_imine.py`/`_nitrile.py` would see for
the plain neutral parent), delegate to that module's own existing naming
logic unchanged, then append the oxide-family suffix word. This project
has no separate acyl-group-name-reuse module (see `_carbenium.py`'s own
acylium docstring for the analogous radical/cation precedent) to lean on
instead.

Explicitly out of scope (raise `UnsupportedStructure`, via whichever
existing module's own validation the stripped molecule falls through to):
anything `_imine.py`/`_nitrile.py` themselves don't support for the
stripped parent (branched imine nitrogen substituent, a ring, more than
one nitrile group, etc.), and this step's own out-of-scope dipole
subtypes: nitrile imide, amine imide (P-74.2.2.2.1.1/P-74.2.1.3 -- the
zwitterionic hydrazinium-ide naming method instead of simple functional-
class 'oxide' naming), and phosphine imides.
"""

from rdkit import Chem

from ._imine import name_imine
from ._nitrile import name_nitrile

# P-66.5.1.2.1: HCN keeps its retained PIN 'formonitrile' (no carbon chain to name).
_NITRILE_OXIDE_RETAINED_OVERRIDES = {"methanenitrile": "formonitrile"}


def _dipole_nitrogen_oxygen(mol):
    """(nitrogen_atom, oxygen_atom) for the sole +1-charged nitrogen
    bonded to a sole -1-charged, degree-1 oxygen, or None if `mol` doesn't
    have exactly this charge pattern (regardless of the rest of its
    shape -- the nitrone-vs-nitrile-oxide distinction is made by the
    caller, from the nitrogen's other bonds)."""
    charged_nitrogens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7 and a.GetFormalCharge() == 1]
    if len(charged_nitrogens) != 1:
        return None
    (nitrogen,) = charged_nitrogens
    charged_oxygens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1]
    if len(charged_oxygens) != 1:
        return None
    (oxygen,) = charged_oxygens
    if oxygen.GetDegree() != 1 or oxygen.GetIsotope() != 0:
        return None
    bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), oxygen.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return None
    return nitrogen, oxygen


def _strip_dipole_oxygen(mol, nitrogen_idx, oxygen_idx):
    """A copy of `mol` with `oxygen_idx` removed and `nitrogen_idx`'s
    charge reset to neutral -- exactly the molecule the plain (non-dipole)
    parent imine/nitrile would be."""
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(oxygen_idx)
    new_nitrogen_idx = nitrogen_idx if nitrogen_idx < oxygen_idx else nitrogen_idx - 1
    rw.GetAtomWithIdx(new_nitrogen_idx).SetFormalCharge(0)
    rw.GetAtomWithIdx(new_nitrogen_idx).SetNoImplicit(False)
    new_mol = rw.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol


def has_nitrone_shape(mol) -> bool:
    """True if `mol` has the dipole N+/O- charge pattern and the charged
    nitrogen also carries a C=N double bond (imine-shaped once the dipole
    oxygen is set aside) -- a nitrone."""
    pair = _dipole_nitrogen_oxygen(mol)
    if pair is None:
        return False
    nitrogen, _ = pair
    return any(
        bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(nitrogen).GetAtomicNum() == 6
        for bond in nitrogen.GetBonds()
    )


def has_nitrile_oxide_shape(mol) -> bool:
    """True if `mol` has the dipole N+/O- charge pattern and the charged
    nitrogen also carries a C#N triple bond (nitrile-shaped once the
    dipole oxygen is set aside) -- a nitrile oxide."""
    pair = _dipole_nitrogen_oxygen(mol)
    if pair is None:
        return False
    nitrogen, _ = pair
    return any(
        bond.GetBondTypeAsDouble() == 3.0 and bond.GetOtherAtom(nitrogen).GetAtomicNum() == 6
        for bond in nitrogen.GetBonds()
    )


def name_nitrone(mol) -> str:
    nitrogen, oxygen = _dipole_nitrogen_oxygen(mol)
    stripped = _strip_dipole_oxygen(mol, nitrogen.GetIdx(), oxygen.GetIdx())
    return f"{name_imine(stripped)} N-oxide"


def _is_bare_hydrogen_nitrile(mol) -> bool:
    """True if `mol` is exactly HC#N (fulminic acid's stripped parent) --
    `_nitrile.py` itself treats a zero-carbon-neighbor nitrile carbon as
    out of scope (no substitutive '-nitrile' suffix name exists for a bare
    terminus), so this shape is special-cased here directly rather than
    reached through that module (P-61.10's own retained 'formonitrile'
    name for this exact parent)."""
    if mol.GetNumAtoms() != 2:
        return False
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    if len(carbons) != 1:
        return False
    (carbon,) = carbons
    (bond,) = carbon.GetBonds()
    return bond.GetBondTypeAsDouble() == 3.0 and bond.GetOtherAtom(carbon).GetAtomicNum() == 7


def name_nitrile_oxide(mol) -> str:
    nitrogen, oxygen = _dipole_nitrogen_oxygen(mol)
    stripped = _strip_dipole_oxygen(mol, nitrogen.GetIdx(), oxygen.GetIdx())
    if _is_bare_hydrogen_nitrile(stripped):
        return "formonitrile oxide"
    parent_name = name_nitrile(stripped)
    parent_name = _NITRILE_OXIDE_RETAINED_OVERRIDES.get(parent_name, parent_name)
    return f"{parent_name} oxide"
