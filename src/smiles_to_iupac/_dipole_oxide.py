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

import re

from rdkit import Chem

from ._common import UnsupportedStructure
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


_CHALCOGEN_TERM = {8: "oxide", 16: "sulfide", 34: "selenide", 52: "telluride"}


def _nitrile_oxide_groups(mol):
    """[(nitrogen, chalcogen)] of every R-C#N(+)-X(-) or R-C#N(=X) group (P-66.5.4.1)."""
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2:
            continue
        triple = chalcogen = False
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            order = bond.GetBondTypeAsDouble()
            if other.GetAtomicNum() == 6 and order == 3.0:
                triple = True
            elif other.GetAtomicNum() in _CHALCOGEN_TERM and other.GetDegree() == 1 and not other.GetIsotope():
                charged = atom.GetFormalCharge() == 1 and other.GetFormalCharge() == -1 and order == 1.0
                neutral = not atom.GetFormalCharge() and not other.GetFormalCharge() and order == 2.0
                if charged or neutral:
                    chalcogen = other.GetIdx()
        if triple and chalcogen is not False:
            found.append((atom.GetIdx(), chalcogen))
    return found


def has_nitrile_oxide_shape(mol) -> bool:
    """True if every charge of `mol` belongs to a nitrile oxide or chalcogen analogue (a salt names its anion first)."""
    groups = _nitrile_oxide_groups(mol)
    if not groups or len(Chem.GetMolFrags(mol)) != 1:
        return False
    own = {atom for pair in groups for atom in pair}
    return not any(a.GetFormalCharge() and a.GetIdx() not in own for a in mol.GetAtoms())


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
    """P-66.5.4.1 method (1): the name of the nitrile with the term 'oxide', 'sulfide', 'selenide' or 'telluride'; the
    nitrile is the parent although a zwitterion outranks esters, amides and acids only in the class order."""
    from ._common import multiplied_word
    from ._polyfunctional import FORCED_PRINCIPAL, name_polyfunctional

    groups = _nitrile_oxide_groups(mol)
    terms = {_CHALCOGEN_TERM[mol.GetAtomWithIdx(x).GetAtomicNum()] for _, x in groups}
    if len(terms) != 1:
        raise UnsupportedStructure("nitrile oxides of different chalcogens are not supported yet")
    stripped = Chem.RWMol(mol)
    for n, _ in groups:
        stripped.GetAtomWithIdx(n).SetFormalCharge(0)
        stripped.GetAtomWithIdx(n).SetNoImplicit(False)
    for x in sorted((x for _, x in groups), reverse=True):
        stripped.RemoveAtom(x)
    stripped = stripped.GetMol()
    Chem.SanitizeMol(stripped)
    if _is_bare_hydrogen_nitrile(stripped):
        parent_name = "formonitrile"
    else:
        token = FORCED_PRINCIPAL.set("nitrile")
        try:
            parent_name = name_polyfunctional(stripped)
        except UnsupportedStructure:
            parent_name = name_nitrile(stripped)
        finally:
            FORCED_PRINCIPAL.reset(token)
        parent_name = _NITRILE_OXIDE_RETAINED_OVERRIDES.get(parent_name, parent_name)
    return f"{parent_name} {multiplied_word(len(groups), terms.pop())}"


_YLIDENE = {8: "oxo", 16: "sulfanylidene", 34: "selanylidene", 52: "tellanylidene"}
_IODO_ONLY_PREFIX = re.compile(r"(?P<head>(?:[a-z]+ )?)(?P<locant>\d+(?:,\d+)*-)?iodo(?P<tail>[a-z].*)")


def has_nitrile_oxide_prefix_shape(mol) -> bool:
    """A nitrile oxide group beside a senior class (an anion or a salt): it is cited as a prefix (P-66.5.4.2)."""
    groups = _nitrile_oxide_groups(mol)
    return len(groups) == 1 and not has_nitrile_oxide_shape(mol)


def name_nitrile_oxide_prefix(mol) -> str:
    """'4-[(oxo-lambda5-azanylidyne)methyl]benzoate': the group is replaced by iodine, the structure named, and the
    iodo prefix, the only one, exchanged for the prefix of the group."""
    (nitrogen, chalcogen), = _nitrile_oxide_groups(mol)
    carbon = next(n for n in mol.GetAtomWithIdx(nitrogen).GetNeighbors() if n.GetAtomicNum() == 6)
    hosts = [n for n in carbon.GetNeighbors() if n.GetIdx() != nitrogen]
    if len(hosts) != 1:
        raise UnsupportedStructure("a nitrile oxide group without a single attachment is not supported as a prefix")
    probe = Chem.RWMol(mol)
    atom = probe.GetAtomWithIdx(carbon.GetIdx())
    atom.SetAtomicNum(53)
    atom.SetNoImplicit(False)
    for idx in sorted((nitrogen, chalcogen), reverse=True):
        probe.RemoveAtom(idx)
    probe = probe.GetMol()
    Chem.SanitizeMol(probe)
    from .core import smiles_to_iupac

    name = smiles_to_iupac(Chem.MolToSmiles(probe))
    match = _IODO_ONLY_PREFIX.fullmatch(name)
    if match is None:
        raise UnsupportedStructure("a nitrile oxide prefix beside other substituents is not supported yet")
    ylidene = _YLIDENE[mol.GetAtomWithIdx(chalcogen).GetAtomicNum()]
    prefix = f"[({ylidene}-\u03bb5-azanylidyne)methyl]"
    return f"{match.group('head')}{match.group('locant') or ''}{prefix}{match.group('tail')}"
