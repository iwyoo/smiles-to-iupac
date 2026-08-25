"""Naming of [2.2]paracyclophane by hardcoded retained-name recognition of
the unsubstituted parent hydride only, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-26 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): phane
  nomenclature names a ring assembly of a macrocycle ("...phane") built
  from component rings ("superatoms") joined by bridges, using composite
  locants and multiplying prefixes on the component-ring names (e.g.
  'benzena' for a benzene component). The commonly used bracket name
  "[2.2]paracyclophane" is only a retained/informal name; the actual PIN,
  confirmed directly against the Blue Book text (P-26.3.2.1 and P-26.4.1.4
  both give the identical string, each explicitly marked "(PIN)"), is
  "1,4(1,4)-dibenzenacyclohexaphane": simplifying each benzene ring to one
  superatom plus the four remaining -CH2- atoms gives a 6-membered
  macrocycle ('cyclohexaphane'); both superatoms sit at 1,4 of that
  macrocycle and are themselves attached to the macrocycle at their own
  1,4 positions (the composite locant '(1,4)'); two identical benzene
  components give the 'di-' multiplying prefix + 'benzena' (P-26.2.3.1 /
  P-26.2.2.1).
- Since this module's only job is recognizing the exact unsubstituted
  parent (no locants to assign, no macrocycle-size or composite-locant
  algorithm needed), an exact whole-molecule match against the canonical
  structure is both sufficient and simplest, mirroring
  `_peri_fused_aromatic.py`/`_branched_fused_aromatic.py`/
  `_heteroaromatic_fused.py`/`_fullerene.py`: any substituent, different
  bridge length, or different substitution pattern (metacyclophane, a
  single-bridge phane, etc.) changes the canonical SMILES and so is
  correctly left unmatched, falling through to `core.py`'s other dispatch
  branches (which all currently raise `UnsupportedStructure` for these
  shapes, since no general phane algorithm exists yet).

Formula cross-checked: C16H16 (PubChem CID 74210; Wikipedia's
"(2.2)Paracyclophane" article independently gives the same PIN). PubChem's
own *computed* "IUPACName" property for this compound is a von Baeyer
bridged-ring name, not a phane name -- it doesn't implement P-26 phane
nomenclature at all, so it isn't usable for verifying the name this module
returns; only the Blue Book's own worked examples are.

Explicitly out of scope: any substituted derivative, any other bridge
length/count (`[3.3]paracyclophane`, single-bridge `[n]paracyclophane`
etc.), `metacyclophane`/other substitution patterns, and any general phane
naming algorithm (composite locants, macrocycle-size calculation) --
`has_paracyclophane_name` simply returns False for all of these, so
`core.py`'s existing dispatch handles them (and continues to raise
`UnsupportedStructure`, unchanged).
"""

from rdkit import Chem

_NAME = "1,4(1,4)-dibenzenacyclohexaphane"
_CANONICAL_SMILES = Chem.CanonSmiles("C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2")


def has_paracyclophane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _CANONICAL_SMILES


def name_paracyclophane(mol) -> str:
    return _NAME
