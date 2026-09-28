"""Naming the best-known *peri*-fused mancude (aromatic) ring systems --
pyrene, acenaphthylene, fluoranthene, aceanthrylene, acephenanthrylene,
coronene, and perylene -- by hardcoded retained-name recognition of the
unsubstituted parent hydride only, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-25.1.2 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  each of these is a retained name (P-25.3.1.3 defines "peri-fusion" -- a
  ring atom shared by three or more rings -- as distinct from the
  ortho-fusion `_aromatic.py` handles; every compound here has such an
  atom). Fluoranthene/aceanthrylene/acephenanthrylene are each exactly a
  plain benzo ring ortho-fused onto acenaphthylene at one of its three
  distinct eligible peripheral bonds (confirmed by building each
  candidate from `acenaphthylene`'s own reference SMILES and comparing
  InChI against the PubChem-retrieved structure for all three names --
  P-25.3.1.3's general fusion-locant-letter mechanism `_pyrene_fusion.py`
  etc. use for other bases can't name these, precisely because a retained
  name always outranks the systematic "benzo[x]acenaphthylene"
  construction once one exists for that exact shape, the same reason
  `_anthracene_fusion.py`/`_phenanthrene_fusion.py` exclude the letters
  that collide with tetracene/chrysene).
- Deliberately not built on `_aromatic.py`'s general ortho-fused-chain
  machinery (periphery walk, orientation, numbering): that algorithm
  assumes no ring atom is shared by 3+ rings (`find_aromatic_fused_core`'s
  own docstring), which is false for every compound here by definition,
  and acenaphthylene's five-membered ring (a non-aromatic -CH=CH- bridge
  across naphthalene's 1,8 positions, confirmed via RDKit's own
  aromaticity perception -- its two bridge carbons are non-aromatic even
  though the whole system is a mancude ring assembly per P-25.1.2) doesn't
  even satisfy `find_aromatic_fused_core`'s "every SSSR ring is a
  6-membered aromatic ring" precondition, so it would never reach that
  module's dispatch at all. The same non-aromatic-bridge shape recurs,
  unchanged, in fluoranthene/aceanthrylene/acephenanthrylene's own
  five-membered ring.
- Since this module's only job is recognizing the exact unsubstituted
  parent (no locants to assign, no numbering needed), an exact whole-
  molecule match against each retained name's canonical structure is both
  sufficient and simplest: any substituent, added/removed hydrogen, or
  additional fused ring changes the canonical SMILES and so is correctly
  left unmatched, falling through to `core.py`'s other (already-correct)
  dispatch branches -- which currently all end up raising
  `UnsupportedStructure` for these shapes, since neither substituted nor
  combined peri-fused systems are supported yet.

Formulas cross-checked: pyrene C16H10, acenaphthylene C12H8, fluoranthene/
aceanthrylene/acephenanthrylene each C16H10 (PubChem CID 9154/107781/9143).

coronene and perylene are two more P-25.1.1 PINs of this same "cyclic
peri-fusion graph, no fixed-numbering sub-piece to reuse" shape
(`tmp/bluebook/P2.txt` lines 2956/2962 list both as PINs) -- confirmed
against PubChem CID 9115 (coronene, C24H12) and CID 9142 (perylene,
C20H12). Both were previously unsupported (`UnsupportedStructure`),
raised by the same "peri-fused... not supported" guard this module's
dict already resolves for the other five names.

Explicitly out of scope: any substituted derivative, any peri-fused ring
system other than these seven exact compounds, and combining any of them
with another fused/bridged/spiro system. `has_retained_peri_fused_name`
simply returns False for all of these, so `core.py`'s existing dispatch
handles them (and continues to raise `UnsupportedStructure`, unchanged).
"""

from rdkit import Chem

_RETAINED_NAME_SMILES = {
    "pyrene": "c1cc2ccc3cccc4ccc(c1)c2c34",
    "acenaphthylene": "C1=Cc2cccc3cccc1c23",
    "fluoranthene": "C1=CC=C2C(=C1)C3=CC=CC4=C3C2=CC=C4",
    "aceanthrylene": "C1=CC=C2C3=C4C(=CC=CC4=CC2=C1)C=C3",
    "acephenanthrylene": "C1=CC=C2C(=C1)C=C3C=CC4=C3C2=CC=C4",
    "coronene": "C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7",
    "perylene": "C1=CC2=C3C(=C1)C4=CC=CC5=C4C(=CC=C5)C3=CC=C2",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for name, smiles in _RETAINED_NAME_SMILES.items()}


def has_retained_peri_fused_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_retained_peri_fused(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
