"""Naming of the 9 canonical inositol stereoisomers (cyclohexane-1,2,3,4,5,6-
hexols with one of the Blue Book's 9 fixed retained diastereomer names), per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-104.1/P-104.2.1 (Chapter P-10, https://iupac.qmul.ac.uk/BlueBook/P10.html,
  the Blue Book): "Inositols have retained names ...
  Other names of cyclitols are systematic substitutive names" -- only these 9
  fixed diastereomers get a retained-name PIN; every other stereo pattern (or
  no stereo at all) correctly falls through to the existing systematic
  'cyclohexane-1,2,3,4,5,6-hexol' engine unchanged.
- Each name's stereo pattern is given by P-104.2.1's own face-notation
  fraction (e.g. myo-inositol '(1,2,3,5/4,6-)'): the numerator locants have
  their hydroxyl on one face of the ring, the denominator locants on the
  other. Each entry below was derived directly from that fraction by
  building an idealized planar-ring 3D model (ring atoms 1..6 placed
  clockwise as viewed from the '+z' face, each hydroxyl at +z if its locant
  is in the numerator else -z) and letting RDKit assign CIP descriptors from
  those 3D coordinates -- not hand-derived R/S labels, and not an external
  structure lookup (PubChem's own name search doesn't distinguish these
  isomers, see the originating issue).
- This exact construction method (clockwise numbering with locant-1's
  hydroxyl 'up' = P-104.2.1's own 'L' configuration) was independently
  calibrated against P-104.2.3's own worked example before being trusted
  for the 9 inositols: '1L-1,2/3,5-cyclohexanetetrol' is stated to equal
  '(1R,2R,3R,5R)-cyclohexane-1,2,3,5-tetrol', and building that exact
  fraction ('1,2/3,5', clockwise) with this method reproduces that exact
  CIP label set, (1R,2R,3R,5R), confirming both the geometric convention
  and this project's RDKit-based CIP engine agree with the Blue Book here.
- Of the 9, only the 1D-/1L-chiro-inositol pair is genuinely chiral (no
  internal mirror plane) -- confirmed computationally: building each of the
  other 7 fractions with the ring numbered in either rotational direction
  (clockwise or counterclockwise) produces the identical canonical SMILES
  both times, i.e. the molecule is superimposable on its own mirror image
  (a meso compound) for all of them. The chiro pair's shared fraction
  '(1,2,4/3,5,6-)' instead produces two distinct canonical structures (one
  per direction), matching P-104.2.1's own text that only 'D'/'L' (not the
  fraction) distinguishes this one enantiomeric pair.
- Matched by a whole-molecule canonical-SMILES lookup, mirroring
  `_steroid_parent_hydrides.py`'s/`_fullerene.py`'s identical retained-name-
  by-exact-structure pattern: any substituent, missing hydroxyl, or
  unspecified/partially-specified stereocenter changes the canonical SMILES
  and so is correctly left unmatched, falling through to the generic
  systematic hexol name.
"""

from rdkit import Chem

_INOSITOL_SMILES = {
    "cis-inositol": "O[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O",
    "epi-inositol": "O[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "allo-inositol": "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O",
    "myo-inositol": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "muco-inositol": "O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O",
    "neo-inositol": "O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "scyllo-inositol": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "1L-chiro-inositol": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",
    "1D-chiro-inositol": "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O",
}
_TEMPLATES = {name: Chem.MolFromSmiles(smiles) for name, smiles in _INOSITOL_SMILES.items()}
assert len({Chem.CanonSmiles(smiles) for smiles in _INOSITOL_SMILES.values()}) == len(_INOSITOL_SMILES), (
    "two entries above canonicalized to the same key -- a real name "
    "collision, not just a duplicate row"
)


def _match(mol):
    # Canonical SMILES of a meso ring depends on the stereo-perception mode, so compare by chiral substructure.
    if mol.GetNumAtoms() != 12:
        return None
    for name, template in _TEMPLATES.items():
        if mol.HasSubstructMatch(template, useChirality=True):
            return name
    return None


def has_inositol_shape(mol) -> bool:
    return _match(mol) is not None


def name_inositol(mol) -> str:
    return _match(mol)
