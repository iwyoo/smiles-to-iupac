"""Naming of triphenylene -- the best-known *branched* ortho-fused mancude
(aromatic) ring system, where one ring is fused to three others rather than
at most two -- by hardcoded retained-name recognition of the unsubstituted
parent hydride only, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-25.1.2 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  'triphenylene' is a retained name.
- This is a different axis from `_peri_fused_aromatic.py`'s pyrene/
  acenaphthylene: those have a single *atom* shared by three or more
  rings (P-25.3.1.3's peri-fusion); triphenylene's four rings are all
  plain ortho-fused (each fusion is a shared *bond*, two atoms), but its
  ring-*adjacency graph* branches -- the central ring is fused to all
  three others, none of which touch each other -- rather than forming the
  simple chain `_aromatic.py`'s `_ring_path_order` requires (that
  function explicitly raises `UnsupportedStructure` naming triphenylene
  as its own branched-topology example). So, like the peri-fused module,
  this bypasses `_aromatic.py`'s chain-numbering machinery entirely rather
  than trying to extend it.
- Since this module's only job is recognizing the exact unsubstituted
  parent (no locants to assign, no numbering needed), an exact whole-
  molecule match against triphenylene's canonical structure is both
  sufficient and simplest -- mirrors `_peri_fused_aromatic.py`'s approach
  exactly; any substituent or different ring-fusion shape changes the
  canonical SMILES and so is correctly left unmatched, falling through to
  `core.py`'s other (already-correct) dispatch branches.

Formula cross-checked: triphenylene C18H12; ring-adjacency degree
sequence [1, 1, 1, 3] (one ring fused to the other three, which are not
fused to each other) confirmed via RDKit -- the shape `_aromatic.py`'s own
`_ring_path_order` names as its motivating branched-topology example.

Explicitly out of scope: any substituted derivative, and any other
branched-fusion ring system (general benzo[x,y-z] combinations with a
branch point are `general-fusion-naming.md`'s territory, not this
module's). `has_retained_branched_fused_name` returns False for both, so
`core.py`'s existing dispatch handles them (and continues to raise
`UnsupportedStructure`, unchanged).
"""

from rdkit import Chem

_RETAINED_NAME_SMILES = {
    "triphenylene": "c1ccc2c(c1)c1ccccc1c1ccccc21",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for name, smiles in _RETAINED_NAME_SMILES.items()}


def has_retained_branched_fused_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_retained_branched_fused(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
