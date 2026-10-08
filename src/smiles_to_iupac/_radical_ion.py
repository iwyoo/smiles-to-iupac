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
  worked example `benzenaminiumyl (PIN)`, the Blue Book
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
  `propyloxidaniumyl (PIN)`, the Blue Book -- the exact
  same "protonate to the full-valence cation, then remove one radical H"
  pattern as `aminiumyl` above, just for O/S instead of N. Confirmed this
  session:

  - `C[OH+]`  (radical=1) -> reconstructs to `C[OH2+]`  -> "methyloxidanium"
    -> "methyloxidaniumyl".
  - `CC[SH+]` (radical=1) -> reconstructs to `CC[SH2+]` -> "ethylsulfanium"
    (note: `_sulfonium.py`'s own established output is "sulfanium", not
    "sulfonium") -> "ethylsulfaniumyl".

- The `aminyliumyl`/`iminyliumyl`/`amidyliumyl` family (P-73.2.3.2's
  'ylium' cation plus P-75.3.1's radical 'yl', stacked): unlike the two
  families above, the intermediate 'ylium' cation itself
  (e.g. `acetamidylium`, `CC(=O)[NH+]`) is *also* inherently
  radical-carrying under RDKit's valence model (confirmed: `radical=2`,
  not a plain closed-shell cation) -- so a single-H healing step isn't
  enough. Instead this reconstructs all the way back to the *neutral*
  amine/imine/amide (zero the charge, restore enough hydrogens to reach
  the neutral atom's own ordinary valence -- mirrors
  `_characteristic_group_radical_name`'s own delegation to
  `_imine.py`/`_amide.py`/`_amine.py`), then applies the combined string
  transform directly on that neutral name: strip the trailing `e`,
  append `"ylium"`, then append `"yl"`. Confirmed worked example
  `acetamidyliumyl (PIN)`, the Blue Book. Confirmed
  this session:

  - `C[N+]`        (radical=3) -> reconstructs to `CN`       -> "methanamine" -> "methanaminyliumyl".
  - `CC=[N+]`       (radical=2) -> reconstructs to `CC=N`     -> "ethanimine"  -> "ethaniminyliumyl".
  - `CC(=O)[N+]`    (radical=3) -> reconstructs to `CC(=O)N`  -> "ethanamide"  -> "ethanamidyliumyl".

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one radical electron on the charged atom, a coexisting charge
  or radical elsewhere in the molecule, or an isotopically modified atom.
- Any shape where the reconstructed neutral-radical-healed molecule
  doesn't match the corresponding module's own `has_*_shape`/dispatch
  (e.g. a substitution pattern that module itself doesn't support
  standalone).
- Selenium (a `propylselaniumyl`-style analogue, if selenonium is
  separately supported) -- not investigated this session, don't assume
  reachability from the oxygen/sulfur pattern alone.
- P-73.2.3.2's own polyamine/polyimine/polyamide multiplicative cations
  and a divalent '-ylidenylium'-shaped version of the `ylium`+`yl` family
  -- mirrors `_radical.py`'s own P-71.3.3 exclusions.
"""

import re

from rdkit import Chem

from ._amide import has_amide_shape, name_amide
from ._amine import name_amine
from ._ammonium import has_ammonium_shape, name_ammonium
from ._common import UnsupportedStructure
from ._imine import has_simple_imine_shape, name_imine
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
    try:
        name = name_fn(healed)
    except UnsupportedStructure:
        return None
    if name.endswith("anilinium"):
        name = name[: -len("anilinium")] + "benzenaminium"
    return name + "yl"


def _neutral_healed(mol, radical):
    """Reconstruct the fully neutral parent by zeroing the charge and
    restoring enough hydrogens to reach the atom's own ordinary (neutral)
    valence -- a *double* healing step, unlike `_healed`'s single
    radical-electron restoration, since the intermediate 'ylium' cation
    is itself still radical-carrying."""
    idx = radical.GetIdx()
    bond_order_sum = sum(
        mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() for n in radical.GetNeighbors()
    )
    target_valence = {7: 3, 8: 2, 16: 2}.get(radical.GetAtomicNum())
    if target_valence is None:
        return None
    new_h = target_valence - bond_order_sum
    if new_h < 0 or new_h != int(new_h):
        return None

    rw = Chem.RWMol(mol)
    atom = rw.GetAtomWithIdx(idx)
    atom.SetFormalCharge(0)
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(int(new_h))
    atom.SetNumRadicalElectrons(0)
    neutral = rw.GetMol()
    try:
        Chem.SanitizeMol(neutral)
    except (Chem.rdchem.AtomValenceException, Chem.rdchem.KekulizeException):
        return None
    return neutral


def _ylium_yl_radical(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons() != 0]
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() not in (2, 3):
        return None
    (radical,) = radicals
    if radical.GetAtomicNum() != 7 or radical.GetFormalCharge() != 1 or radical.GetIsotope() != 0:
        return None
    if radical.GetIsAromatic():
        return None
    if any(
        atom.GetIdx() != radical.GetIdx() and (atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons() != 0)
        for atom in mol.GetAtoms()
    ):
        return None

    neutral = _neutral_healed(mol, radical)
    if neutral is None:
        return None

    try:
        if has_simple_imine_shape(neutral):
            neutral_name = name_imine(neutral)
        elif has_amide_shape(neutral):
            neutral_name = name_amide(neutral)
        else:
            neutral_name = name_amine(neutral)
    except UnsupportedStructure:
        from .core import smiles_to_iupac

        try:
            neutral_name = smiles_to_iupac(Chem.MolToSmiles(neutral))
        except UnsupportedStructure:
            return None
    if not neutral_name.endswith(("amide", "amine", "imine", "aniline")):
        return None
    # P-73.2.3.2: a nitrogen with two missing valences is the hydride-loss cation (acetamidylium); the third is a radical
    return neutral_name[:-1] + "ylium" + ("yl" if radical.GetNumRadicalElectrons() == 3 else "")


_RADICAL_SUFFIX = {1: "yl", 2: "ylidene", 3: "ylidyne"}


def _healed_ion_radical(mol):
    """The ion with the radical electrons of its single charged radical atom filled by hydrogens is named, and the
    ionic ending takes the radical suffix: azanide -> azanidyl, hydrazin-1-ide -> hydrazin-1-id-1-yl, aminium -> aminiumyl."""
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if len(radicals) != 1:
        return None
    (radical,) = radicals
    count = radical.GetNumRadicalElectrons()
    charge = radical.GetFormalCharge()
    if abs(charge) != 1 or radical.GetIsotope() or count not in _RADICAL_SUFFIX:
        return None
    if radical.GetAtomicNum() == 7 and charge == 1 and radical.IsInRing() and radical.GetDegree() == 2 and not radical.GetTotalNumHs():
        return None
    editable = Chem.RWMol(mol)
    atom = editable.GetAtomWithIdx(radical.GetIdx())
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(radical.GetTotalNumHs() + count)
    atom.SetNumRadicalElectrons(0)
    healed = editable.GetMol()
    try:
        Chem.SanitizeMol(healed)
        from .core import smiles_to_iupac

        name = smiles_to_iupac(Chem.MolToSmiles(healed))
    except (Chem.rdchem.AtomValenceException, Chem.rdchem.KekulizeException, UnsupportedStructure):
        return None
    ending = "ide" if charge < 0 else "ium"
    if charge > 0 and name.endswith("ide"):
        cation = re.search(r"-(\d+)-ium-\d+(?:,\d+)*-ide$", name)
        return f"{name[:-1]}-{cation.group(1)}-{_RADICAL_SUFFIX[count]}" if cation else None
    if not name.endswith(ending):
        return None
    suffix = _RADICAL_SUFFIX[count]
    lambda_form = _lambda_ide_radical(radical, name, charge, count)
    if lambda_form is not None:
        return lambda_form
    located = re.search(r"-(\d+)-" + ending + "$", name)
    stem = name[:-1] if charge < 0 else name
    if located:
        return f"{stem}-{located.group(1)}{_added_hydrogen(mol, radical, name)}-{suffix}"
    return stem + suffix


def _lambda_ide_radical(radical, name, charge, count):
    """P-75.2.1: a ring chalcogen carrying both the anionic and the radical centre is cited as a lambda-4 parent with
    the 'ide' and 'yl' suffixes: '1λ4-thiiran-1-id-1-yl'."""
    if radical.GetAtomicNum() not in (16, 34, 52) or charge != -1 or count != 1 or not radical.IsInRing():
        return None
    match = re.fullmatch(r"([a-z]+)-(\d+)-uide", name)
    if match is None:
        return None
    return f"{match.group(2)}\u03bb4-{match.group(1)}-{match.group(2)}-id-{match.group(2)}-yl"


def _added_hydrogen(mol, radical, name):
    """P-75.2.4: a ring nitrogen centre beside a ring carbonyl carbon cites the added hydrogen of that oxo group:
    '1-ethyl-2-oxopyridin-1-ium-1(2H)-yl'."""
    oxo = re.findall(r"(?<![\d,])(\d+)-oxo", name)
    if radical.GetAtomicNum() != 7 or not radical.IsInRing() or len(oxo) != 1:
        return ""
    beside = [
        n
        for n in radical.GetNeighbors()
        if n.IsInRing()
        and any(
            m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 2.0
            for m in n.GetNeighbors()
        )
    ]
    return f"({oxo[0]}H)" if len(beside) == 1 else ""


_GROUP_13 = {5, 13, 31, 49, 81}


def _ammonium_group(mol, nitrogen, centre):
    """'N,N-diethylethanaminiumyl': the substituent group of a quaternary nitrogen, named as the ammonium cation that
    has a hydrogen in place of the bond to `centre`, with the radical suffix."""
    keep, stack = {nitrogen}, [nitrogen]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetIdx() != centre and n.GetIdx() not in keep:
                keep.add(n.GetIdx())
                stack.append(n.GetIdx())
    editable = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(index)
    new_index = sorted(keep).index(nitrogen)
    target = editable.GetAtomWithIdx(new_index)
    target.SetNoImplicit(True)
    target.SetNumExplicitHs(1)
    cation = editable.GetMol()
    Chem.SanitizeMol(cation)
    if not has_ammonium_shape(cation):
        return None
    name = name_ammonium(cation)
    if name.endswith("anilinium"):
        name = name[: -len("anilinium")] + "benzenaminium"
    return name + "yl"


def _group13_zwitterion(mol):
    """A group 13 centre carrying the radical and the negative charge, with ammonium groups beside organyl groups
    (P-75.2.3, P-75.4): '(N,N-diethylethanaminiumyl)boranuidyl'."""
    from ._common import adjacency, halogen_substituents
    from ._hetero_prefixes import MONONUCLEAR_HYDRIDES
    from ._substituents import format_mononuclear_prefixes, name_branch

    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if len(radicals) != 1 or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetIsotope() for a in mol.GetAtoms()):
        return None
    (centre,) = radicals
    if centre.GetAtomicNum() not in _GROUP_13 or centre.GetFormalCharge() != -1 or centre.GetNumRadicalElectrons() != 1 or centre.IsInRing():
        return None
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge() > 0]
    if not cations or any(a.GetAtomicNum() != 7 or a.GetFormalCharge() != 1 or a.GetIdx() not in {n.GetIdx() for n in centre.GetNeighbors()} for a in cations):
        return None
    if len(cations) + sum(1 for a in mol.GetAtoms() if a.GetFormalCharge() < 0) != 2 * len(cations):
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = []
    for neighbour in centre.GetNeighbors():
        if neighbour.GetFormalCharge() > 0:
            group = _ammonium_group(mol, neighbour.GetIdx(), centre.GetIdx())
            if group is None:
                return None
            entries.append((f"({group})", False))
        elif neighbour.GetAtomicNum() == 6:
            entries.append(name_branch(graph, neighbour.GetIdx(), centre.GetIdx(), halogens, aromatic, mol=mol))
        else:
            return None
    stem = MONONUCLEAR_HYDRIDES[centre.GetAtomicNum()][0][:-1] + "uidyl"
    return format_mononuclear_prefixes(entries) + stem


def has_radical_ion_shape(mol) -> bool:
    """True if `mol` matches `_ionic_suffix_radical`'s or
    `_ylium_yl_radical`'s own P-75.3.1/.2 shape. Used by `core.py` to
    route here ahead of `has_radical_shape`, whose own broader "any
    nonzero radical electron count" check would otherwise claim this
    charge+radical combination first and misroute it into the
    plain-radical dispatch."""
    return (
        _ionic_suffix_radical(mol) is not None
        or _ylium_yl_radical(mol) is not None
        or _group13_zwitterion(mol) is not None
        or _healed_ion_radical(mol) is not None
    )


def name_radical_ion(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    name = _ionic_suffix_radical(mol)
    if name is not None:
        return name
    name = _ylium_yl_radical(mol) or _group13_zwitterion(mol) or _healed_ion_radical(mol)
    if name is None:
        raise UnsupportedStructure(
            "only the aminiumyl/oxidaniumyl/sulfaniumyl/aminyliumyl/"
            "iminyliumyl/amidyliumyl radical cations are supported yet "
            "(P-75.3.1/.2)"
        )
    return name
