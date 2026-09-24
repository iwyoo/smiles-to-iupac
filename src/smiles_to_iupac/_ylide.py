"""Naming of nitrogen ylides (R3N+-CR2-, P-74.2.1.1.1's Method (1)), per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-74.2.1.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf,
  ~2867-2874): a nitrogen ylide is named as a zwitterion -- the anionic
  carbon is the parent, suffixed '-ide' (P-72.2.2.1's general carbanion
  construction, elision of the parent hydride's final 'e'), and the
  cationic nitrogen group is cited as a substituent prefix, its own name
  derived from the corresponding neutral amine (`_amine.py`'s own
  secondary/tertiary construction, reused directly -- the same machinery
  `_ammonium.py`'s own quaternary-ammonium path already reuses) with
  'amine' swapped for 'aminium' (P-73.1.1.2's hydron-addition rule) plus a
  trailing '-yl' for the point of attachment. Confirmed worked example
  `2-(N,N-dimethylmethanaminiumyl)propan-2-ide (PIN)`.
- Scoped here to a *mononuclear* anion carbon only (the ylide carbon's
  sole bond is the ylide bond itself to the ammonium nitrogen, always
  giving the constant parent name 'methanide') -- the worked example's own
  branched-chain anion case (`propan-2-ide`, the charge sitting at an
  *internal* chain position, cited with its own locant) needs its own
  chain-numbering/locant-elision research this session didn't confirm
  against a second primary-source example, so it's deliberately narrower
  than the milestone's own general "R3N+-CR2-" framing -- a follow-up,
  not guessed at here. Confirmed via a real PubChem structure (CID
  12436722, `[CH2-][N+](C)(C)C`) for exactly this mononuclear shape.
  `_ammonium.py`'s own quaternary-ammonium check already explicitly
  excludes this shape (a second charged atom elsewhere) rather than
  silently mis-naming it -- this module is the follow-up that actually
  names it instead of just rejecting it.

- P-74.2.1.1.2/.3/.4 (phosphorus/oxygen/sulfur ylides, R3P+-CR2-/
  R2O+-CR2-/R2S+-CR2-): the identical zwitterionic method, one heteroatom
  over -- the cationic group's own name is derived from the corresponding
  neutral cation module (`_phosphonium.py`/`_oxonium.py`/`_sulfonium.py`)
  instead of `_amine.py`, by isolating the cation atom and its own
  substituents into a standalone fragment (cutting the ylide bond,
  adding an extra hydrogen -- the same "cut a bond, add an H" trick
  `_zwitterion.py`'s own `_ammonium_prefix` already uses) and naming that
  fragment with `name_phosphonium`/`name_oxonium`/`name_sulfonium`
  directly, then appending '-yl'. Confirmed this session by running each
  namer directly on the isolated fragment: `name_phosphonium` on
  `[P+](C)(C)C` -> `"trimethylphosphanium"`, `name_oxonium` on `[O+](C)C`
  -> `"dimethyloxidanium"`, `name_sulfonium` on `[S+](C)C` ->
  `"dimethylsulfanium"` -- each `+ "yl"` matches the milestone's own
  cited worked examples' cation half exactly. Same mononuclear
  ('methanide') anion-carbon restriction as the nitrogen case above (see
  that shape's own note on the still-open branched-anion follow-up,
  applicable to all four heteroatoms alike, not repeated per-element).
  Oxygen/sulfur cations only ever carry 2 substituents (trivalent once
  charged, one bond being the ylide bond itself), phosphorus/nitrogen up
  to 3 (tetravalent once charged) -- `_MAX_CATION_DEGREE` below encodes
  this per element.

Explicitly out of scope (raise `UnsupportedStructure`):
- A branched/chain carbanion parent (anion carbon bonded to anything
  besides the ylide bond itself) -- see above.
- A cation atom not fully substituted to its own maximum valence (an N-H/
  P-H ylide, or an O/S ylide with only 1 substituent, isn't yet confirmed
  against a primary-source example), any cation substituent not shaped
  like a plain saturated branch already within `_amine.py`/
  `_phosphonium.py`/`_oxonium.py`/`_sulfonium.py`'s own existing scope, a
  ring anywhere in the molecule, more than one anionic/cationic center,
  or any heteroatom other than the ylide's own carbon/cation-atom pair.
- Carbene/nitrene dipolar forms (P-74.2.2.3) -- a separate milestone step.
"""

from rdkit import Chem

from ._amine import _name_acyclic_secondary_tertiary_amine
from ._common import UnsupportedStructure, non_single_bonds
from ._numerals import alkane_name
from ._oxonium import name_oxonium
from ._phosphonium import name_phosphonium
from ._sulfonium import name_sulfonium

_MAX_CATION_DEGREE = {7: 4, 15: 4, 8: 3, 16: 3}
_CATION_NAMERS = {15: name_phosphonium, 8: name_oxonium, 16: name_sulfonium}


def _ylide_core(mol, cation_elements):
    """(anion_atom, cation_atom) if `mol` is an ylide-shaped zwitterion
    (see module docstring) whose cation atom is one of `cation_elements`,
    else None. Doesn't itself check the anion carbon's own degree
    (mononuclear-only, this module's own scope) -- that's the caller's
    job, so `has_*_ylide_shape` still correctly returns False for an
    out-of-scope branched anion instead of routing here and only then
    raising."""
    negatives = [a for a in mol.GetAtoms() if a.GetFormalCharge() == -1]
    positives = [a for a in mol.GetAtoms() if a.GetFormalCharge() == 1]
    if len(negatives) != 1 or len(positives) != 1:
        return None
    (anion,) = negatives
    (cation,) = positives
    if anion.GetAtomicNum() != 6 or cation.GetAtomicNum() not in cation_elements:
        return None
    if anion.GetIsotope() != 0 or cation.GetIsotope() != 0:
        return None
    if anion.GetIsAromatic() or cation.GetIsAromatic():
        return None
    bond = mol.GetBondBetweenAtoms(anion.GetIdx(), cation.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return None
    if mol.GetRingInfo().NumRings() != 0:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
    for atom in mol.GetAtoms():
        if atom.GetIdx() in (anion.GetIdx(), cation.GetIdx()):
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0 or atom.GetIsAromatic():
            return None
        if atom.GetAtomicNum() != 6:
            return None
    if any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()):
        return None
    return anion, cation


def has_nitrogen_ylide_shape(mol) -> bool:
    return _ylide_core(mol, {7}) is not None


def name_nitrogen_ylide(mol) -> str:
    anion, cation = _ylide_core(mol, {7})
    if anion.GetDegree() != 1:
        raise UnsupportedStructure(
            "a nitrogen-ylide anion carbon bonded to anything other than "
            "the ylide nitrogen itself is not supported yet (a branched/"
            "chain carbanion parent)"
        )
    if cation.GetDegree() != 4:
        raise UnsupportedStructure(
            "a nitrogen-ylide ammonium nitrogen not bonded to exactly "
            "three other carbon substituents (plus the ylide bond itself) "
            "is not supported yet"
        )
    n_carbons = tuple(n.GetIdx() for n in cation.GetNeighbors() if n.GetIdx() != anion.GetIdx())
    if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in n_carbons):
        raise UnsupportedStructure("an ammonium nitrogen substituent other than carbon is not supported yet")

    bonds = [b for b in non_single_bonds(mol) if b[2] in (2.0, 3.0)]
    amine_name = _name_acyclic_secondary_tertiary_amine(mol, cation.GetIdx(), n_carbons, bonds)
    aminium_name = amine_name[:-1] + "ium"

    anion_name = alkane_name(1)[:-1] + "ide"
    return f"({aminium_name}yl){anion_name}"


def has_pos_ylide_shape(mol) -> bool:
    return _ylide_core(mol, {15, 8, 16}) is not None


def name_pos_ylide(mol) -> str:
    anion, cation = _ylide_core(mol, {15, 8, 16})
    if anion.GetDegree() != 1:
        raise UnsupportedStructure(
            "a phosphorus/oxygen/sulfur-ylide anion carbon bonded to "
            "anything other than the ylide bond itself is not supported "
            "yet (a branched/chain carbanion parent)"
        )
    max_degree = _MAX_CATION_DEGREE[cation.GetAtomicNum()]
    if cation.GetDegree() != max_degree:
        raise UnsupportedStructure(
            "a phosphorus/oxygen/sulfur-ylide cation atom not fully "
            "substituted to its own maximum valence is not supported yet"
        )
    other_neighbors = [n.GetIdx() for n in cation.GetNeighbors() if n.GetIdx() != anion.GetIdx()]
    if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in other_neighbors):
        raise UnsupportedStructure("a cation substituent other than carbon is not supported yet")

    rw = Chem.RWMol(mol)
    cation_idx = cation.GetIdx()
    rw.RemoveBond(cation_idx, anion.GetIdx())
    rw.GetAtomWithIdx(cation_idx).SetNoImplicit(True)
    rw.GetAtomWithIdx(cation_idx).SetNumExplicitHs(mol.GetAtomWithIdx(cation_idx).GetTotalNumHs() + 1)
    isolated_mol = rw.GetMol()
    Chem.SanitizeMol(isolated_mol)
    mapping = []
    fragments = Chem.GetMolFrags(isolated_mol, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=mapping)
    (cation_fragment,) = (
        fragment for fragment, atom_indices in zip(fragments, mapping) if cation_idx in atom_indices
    )
    cation_name = _CATION_NAMERS[cation.GetAtomicNum()](cation_fragment)

    anion_name = alkane_name(1)[:-1] + "ide"
    return f"({cation_name}yl){anion_name}"
