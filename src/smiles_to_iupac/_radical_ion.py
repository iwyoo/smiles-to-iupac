"""Naming of radical ions on ionic suffix groups (P-75.3.1), per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-75 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf): a
  radical ion is a species with both a radical center and an ionic
  (charged) center, which may sit on the same atom. P-75.3.1 ("Radical
  ions on ionic suffix groups") covers the case where the ionic center is
  already named by an existing ionic suffix group -- the radical suffix
  'yl' is simply added on top, with elision of the final 'e' where the
  ionic parent name ends in one.

- The `aminiumyl` family: an ordinary ammonium cation (`_ammonium.py`'s
  own shape, P-73.1.1.2 -- a nitrogen formed by adding a hydron to a
  neutral amine) with one hydrogen further removed as a radical. Confirmed
  worked example `benzenaminiumyl (PIN)`, `tmp/bluebook/P7.txt` ~3478-3496
  -- built from `benzenaminium` (the ordinary anilinium/ammonium cation)
  plus the radical `-yl` suffix; no elision, since `_ammonium.py`'s own
  names never end in 'e'.

  Discovered while implementing #947/M2 (cation centers on characteristic
  groups): the amine/imine/amide- and hydroxy/chalcogen-derived cations
  P-73.2.3.2/.3/.4 describe (formed by removing a *hydride ion* from an
  already-saturated neutral group) turn out to be structurally
  unreachable as plain closed-shell cations under this project's
  RDKit-based molecular representation -- any real SMILES for that shape
  already carries nonzero radical electrons (RDKit's valence model always
  fills the gap to a charged N/O/S atom's expected valence with radical
  character rather than a lone pair), so it's actually this P-75.3.1
  shape, not a plain P-73.2.3.2 cation.

  Reconstructing the neutral-radical-healed molecule (replacing the
  radical electron with one additional explicit hydrogen, the same
  "strip the marker, sanitize, delegate to the neutral namer" pattern
  `_radical.py`'s own `_characteristic_group_radical_name` and
  `_dipole_oxide.py` already use) turns the radical cation back into an
  ordinary, already-supported `_ammonium.py` shape -- confirmed for
  primary/secondary substitution and an aromatic ring:

  - `C[NH2+]` (radical=1) -> reconstructs to `C[NH3+]` -> "methanaminium"
    -> "methanaminiumyl".
  - `c1ccccc1[NH2+]` (radical=1) -> reconstructs to `c1ccccc1[NH3+]` ->
    "anilinium" -> "aniliniumyl".
  - `C[NH+]C` (radical=1) -> reconstructs to `C[NH2+]C` ->
    "N-methylmethanaminium" -> "N-methylmethanaminiumyl".

- The `oxidaniumyl`/`sulfaniumyl` family (P-75.3.2): the same mechanism
  extended to oxygen/sulfur -- an ordinary oxonium/sulfonium cation
  (`_oxonium.py`/`_sulfonium.py`'s own shape, P-73.1.1.2) with one
  hydrogen further removed as a radical. Confirmed worked example
  `propyloxidaniumyl (PIN)`, `tmp/bluebook/P7.txt` ~3496-3512 -- the exact
  same "protonate to the full-valence cation, then remove one radical H"
  pattern as `aminiumyl` above, just for O/S instead of N. Confirmed this
  session:

  - `C[OH+]`  (radical=1) -> reconstructs to `C[OH2+]`  -> "methyloxidanium"
    -> "methyloxidaniumyl".
  - `CC[SH+]` (radical=1) -> reconstructs to `CC[SH2+]` -> "ethylsulfanium"
    (note: `_sulfonium.py`'s own established output is "sulfanium", not
    "sulfonium") -> "ethylsulfaniumyl".

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one radical electron on the charged atom, a coexisting charge
  or radical elsewhere in the molecule, or an isotopically modified atom.
- Any shape where the reconstructed neutral-radical-healed molecule
  doesn't match the corresponding module's own `has_*_shape` (e.g. a
  substitution pattern that module itself doesn't support standalone).
- The 'ylium'-derived remainder of P-73.2.3.2 (`acetamidyliumyl`-style,
  needing a *double* reconstruction all the way back to the neutral
  amine/imine/amide) -- a separate follow-up step, not this one.
- Selenium (a `propylselaniumyl`-style analogue, if selenonium is
  separately supported) -- not investigated this session, don't assume
  reachability from the oxygen/sulfur pattern alone.
"""

from rdkit import Chem

from ._ammonium import has_ammonium_shape, name_ammonium
from ._common import UnsupportedStructure
from ._oxonium import has_oxonium_shape, name_oxonium
from ._sulfonium import has_sulfonium_shape, name_sulfonium

_IONIC_SUFFIX_FAMILIES = (
    (7, has_ammonium_shape, name_ammonium),
    (8, has_oxonium_shape, name_oxonium),
    (16, has_sulfonium_shape, name_sulfonium),
)


def _healed(mol, radical):
    rw = Chem.RWMol(mol)
    idx = radical.GetIdx()
    atom = rw.GetAtomWithIdx(idx)
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(mol.GetAtomWithIdx(idx).GetTotalNumHs() + 1)
    atom.SetNumRadicalElectrons(0)
    healed = rw.GetMol()
    try:
        Chem.SanitizeMol(healed)
    except (Chem.rdchem.AtomValenceException, Chem.rdchem.KekulizeException):
        return None
    return healed


def _ionic_suffix_radical(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons() != 0]
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() != 1:
        return None
    (radical,) = radicals
    if radical.GetFormalCharge() != 1 or radical.GetIsotope() != 0:
        return None
    family = next((f for f in _IONIC_SUFFIX_FAMILIES if f[0] == radical.GetAtomicNum()), None)
    if family is None:
        return None
    if any(
        atom.GetIdx() != radical.GetIdx() and (atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons() != 0)
        for atom in mol.GetAtoms()
    ):
        return None

    _, has_shape, name_fn = family
    healed = _healed(mol, radical)
    if healed is None or not has_shape(healed):
        return None
    return name_fn(healed) + "yl"


def has_radical_ion_shape(mol) -> bool:
    """True if `mol` matches `_ionic_suffix_radical`'s own P-75.3.1/.2
    shape. Used by `core.py` to route here ahead of `has_radical_shape`,
    whose own broader "any nonzero radical electron count" check would
    otherwise claim this charge+radical combination first and misroute it
    into the plain-radical dispatch."""
    return _ionic_suffix_radical(mol) is not None


def name_radical_ion(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    name = _ionic_suffix_radical(mol)
    if name is None:
        raise UnsupportedStructure(
            "only the aminiumyl/oxidaniumyl/sulfaniumyl radical cations are supported yet (P-75.3.1/.2)"
        )
    return name
