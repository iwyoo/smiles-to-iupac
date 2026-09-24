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

Explicitly out of scope (raise `UnsupportedStructure`):
- A branched/chain carbanion parent (anion carbon bonded to anything
  besides the ylide nitrogen itself) -- see above.
- An ammonium nitrogen with degree other than 4 (an N-H ylide, e.g.
  R2N+(H)-CR2-, isn't yet confirmed against a primary-source example),
  any ammonium substituent not shaped like a plain saturated
  `_amine.py`-supported branch, a ring anywhere in the molecule, more
  than one anionic/cationic center, or any heteroatom other than the
  ylide's own carbon/nitrogen pair.
- Phosphorus/oxygen/sulfur ylides (P-74.2.1.1.2/.3/.4) and carbene/
  nitrene dipolar forms -- each its own separate milestone step.
"""

from rdkit import Chem

from ._amine import _name_acyclic_secondary_tertiary_amine
from ._common import UnsupportedStructure, non_single_bonds
from ._numerals import alkane_name


def _nitrogen_ylide_core(mol):
    """(anion_atom, cation_atom) if `mol` is a nitrogen-ylide-shaped
    zwitterion (see module docstring), else None. Doesn't itself check the
    anion carbon's own degree (mononuclear-only, this module's own
    scope) -- that's `name_nitrogen_ylide`'s job, so `has_nitrogen_ylide_
    shape` still correctly returns False for an out-of-scope branched
    anion instead of routing here and only then raising."""
    negatives = [a for a in mol.GetAtoms() if a.GetFormalCharge() == -1]
    positives = [a for a in mol.GetAtoms() if a.GetFormalCharge() == 1]
    if len(negatives) != 1 or len(positives) != 1:
        return None
    (anion,) = negatives
    (cation,) = positives
    if anion.GetAtomicNum() != 6 or cation.GetAtomicNum() != 7:
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
    return _nitrogen_ylide_core(mol) is not None


def name_nitrogen_ylide(mol) -> str:
    anion, cation = _nitrogen_ylide_core(mol)
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
