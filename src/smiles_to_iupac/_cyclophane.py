"""Naming of [2.2]paracyclophane, [2.2]metacyclophane, and
[1.1.1.1]metacyclophane ("cyclotetrabenzylene") by hardcoded retained-name
recognition of each unsubstituted parent hydride, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-26 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): phane
  nomenclature names a ring assembly of a macrocycle ("...phane") built
  from component rings ("superatoms") joined by bridges, using composite
  locants and multiplying prefixes on the component-ring names (e.g.
  'benzena' for a benzene component). The commonly used bracket names
  "[2.2]paracyclophane"/"[2.2]metacyclophane"/"[1.1.1.1]metacyclophane" are
  only retained/informal names; the actual PINs, confirmed directly against
  the Blue Book text (P-26.3.2.1/P-26.4.1.4 for the first two, P-26.4.1.4's
  own third worked example -- printed immediately alongside the other two,
  all three explicitly marked "(PIN)" -- for the third), are
  "1,4(1,4)-dibenzenacyclohexaphane", "1,4(1,3)-dibenzenacyclohexaphane",
  and "1,3,5,7(1,3)-tetrabenzenacyclooctaphane" respectively: simplifying
  each benzene ring to one superatom plus its remaining -CH2- bridge
  atom(s) gives a 6-membered macrocycle ('cyclohexaphane', two benzena
  superatoms) for the first two, or an 8-membered one ('cyclooctaphane',
  four benzena superatoms) for the third; each superatom's own attachment
  locants on its own benzene ring -- '(1,4)' (para), '(1,3)' (meta, both
  the second and third names) -- are otherwise the only difference between
  the names (P-26.2.3.1/P-26.2.2.1 for the 'di-'/'tetra-'+'benzena'
  multiplying prefix common to all three).
- Since this module's only job is recognizing each exact unsubstituted
  parent (no locants to assign, no macrocycle-size or composite-locant
  algorithm needed), an exact whole-molecule match against each canonical
  structure is both sufficient and simplest, mirroring
  `_peri_fused_aromatic.py`/`_branched_fused_aromatic.py`/
  `_heteroaromatic_fused.py`/`_fullerene.py`: any substituent, different
  bridge length/count, or different substitution pattern changes the
  canonical SMILES and so is correctly left unmatched, falling through to
  `core.py`'s other dispatch branches (which all currently raise
  `UnsupportedStructure` for these shapes, since no general phane algorithm
  exists yet).

Formulas cross-checked. [2.2]Paracyclophane (C16H16): PubChem CID 74210;
Wikipedia's "(2.2)Paracyclophane" article independently gives the same
PIN. [2.2]Metacyclophane (C16H16): PubChem CID 137543 "(2.2)Metacyclophane"
(canonical SMILES independently matches this module's own).
[1.1.1.1]Metacyclophane (C28H24, as of
`tasks/cyclophane-third-case-naming.md`, 2026-08-25): independently built
SMILES resolves via InChIKey lookup to PubChem CID 11740710 (same formula,
same connectivity); also the known parent skeleton of
cyclotetraveratrylene (CID 636089, its octamethoxy derivative), a
well-documented calixarene-precursor macrocycle in the literature. PubChem's
own *computed* "IUPACName" property for any of the three compounds is a von
Baeyer bridged-ring name, not a phane name -- it doesn't implement P-26
phane nomenclature at all, so it isn't usable for verifying any name this
module returns; only the Blue Book's own worked examples are.

Explicitly out of scope: any substituted derivative, any other bridge
length/count/component count, any phane besides these three exact
structures (e.g. `[3.3]paracyclophane`, single-bridge `[n]paracyclophane`),
and any general phane naming algorithm (composite locants, macrocycle-size
calculation) -- `has_cyclophane_name` simply returns False for all of
these, so `core.py`'s existing dispatch handles them (and continues to
raise `UnsupportedStructure`, unchanged).
"""

from rdkit import Chem

_NAMES_BY_SMILES = {
    "C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2": "1,4(1,4)-dibenzenacyclohexaphane",
    "C1CC2=CC(=CC=C2)CCC3=CC=CC1=C3": "1,4(1,3)-dibenzenacyclohexaphane",
    "C1c2cccc(c2)Cc2cccc(c2)Cc2cccc(c2)Cc2cccc1c2": "1,3,5,7(1,3)-tetrabenzenacyclooctaphane",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for smiles, name in _NAMES_BY_SMILES.items()}


def has_cyclophane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_cyclophane(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
