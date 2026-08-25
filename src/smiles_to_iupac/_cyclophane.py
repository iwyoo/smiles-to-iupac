"""Naming of [2.2]paracyclophane and [2.2]metacyclophane by hardcoded
retained-name recognition of each unsubstituted parent hydride, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-26 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): phane
  nomenclature names a ring assembly of a macrocycle ("...phane") built
  from component rings ("superatoms") joined by bridges, using composite
  locants and multiplying prefixes on the component-ring names (e.g.
  'benzena' for a benzene component). The commonly used bracket names
  "[2.2]paracyclophane"/"[2.2]metacyclophane" are only retained/informal
  names; the actual PINs, confirmed directly against the Blue Book text
  (P-26.3.2.1/P-26.4.1.4 for the para isomer, P-26.4.1.4's own paired
  example for the meta isomer -- both explicitly marked "(PIN)"), are
  "1,4(1,4)-dibenzenacyclohexaphane" and "1,4(1,3)-dibenzenacyclohexaphane"
  respectively: simplifying each benzene ring to one superatom plus the
  four remaining -CH2- atoms gives a 6-membered macrocycle
  ('cyclohexaphane'); both superatoms sit at 1,4 of that macrocycle; each
  superatom's own attachment locants on its benzene ring -- '(1,4)' (para)
  or '(1,3)' (meta) -- are the only difference between the two PINs (P-26.2.3.1
  / P-26.2.2.1 for the 'di-'+'benzena' multiplying prefix common to both).
- Since this module's only job is recognizing each exact unsubstituted
  parent (no locants to assign, no macrocycle-size or composite-locant
  algorithm needed), an exact whole-molecule match against each canonical
  structure is both sufficient and simplest, mirroring
  `_peri_fused_aromatic.py`/`_branched_fused_aromatic.py`/
  `_heteroaromatic_fused.py`/`_fullerene.py`: any substituent, different
  bridge length, or different substitution pattern changes the canonical
  SMILES and so is correctly left unmatched, falling through to
  `core.py`'s other dispatch branches (which all currently raise
  `UnsupportedStructure` for these shapes, since no general phane algorithm
  exists yet).

Formulas cross-checked: both C16H16. [2.2]Paracyclophane: PubChem CID
74210; Wikipedia's "(2.2)Paracyclophane" article independently gives the
same PIN. [2.2]Metacyclophane: PubChem CID 137543 "(2.2)Metacyclophane"
(canonical SMILES independently matches this module's own). PubChem's own
*computed* "IUPACName" property for either compound is a von Baeyer
bridged-ring name, not a phane name -- it doesn't implement P-26 phane
nomenclature at all, so it isn't usable for verifying either name this
module returns; only the Blue Book's own worked examples are.

Explicitly out of scope: any substituted derivative, any other bridge
length/count (`[3.3]paracyclophane`, single-bridge `[n]paracyclophane`
etc.), any phane besides these two exact structures, and any general phane
naming algorithm (composite locants, macrocycle-size calculation) --
`has_cyclophane_name` simply returns False for all of these, so
`core.py`'s existing dispatch handles them (and continues to raise
`UnsupportedStructure`, unchanged).
"""

from rdkit import Chem

_NAMES_BY_SMILES = {
    "C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2": "1,4(1,4)-dibenzenacyclohexaphane",
    "C1CC2=CC(=CC=C2)CCC3=CC=CC1=C3": "1,4(1,3)-dibenzenacyclohexaphane",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for smiles, name in _NAMES_BY_SMILES.items()}


def has_cyclophane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_cyclophane(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
