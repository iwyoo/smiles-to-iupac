"""Naming of buckminsterfullerene (C60) by hardcoded recognition of the
exact, unsubstituted cage, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-27 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) and the
  companion fullerene nomenclature document
  (https://iupac.qmul.ac.uk/fullerene/index.html): the icosahedral C60
  cage (12 pentagons + 20 hexagons, every skeletal atom degree 3, the
  "soccer ball"/truncated-icosahedron structure) is named '[60]fullerene'.
- Reference structure taken from PubChem CID 123591
  (buckminsterfullerene)'s own connectivity SMILES, cross-checked here
  via RDKit: 60 all-carbon atoms, every atom degree 3, ring perception of
  exactly 12 five-membered and 20 six-membered rings.
- An exact whole-molecule canonical-SMILES match (mirrors
  `_peri_fused_aromatic.py`'s approach) is used rather than a looser
  shape check (e.g. just the ring-size counts): Euler's formula forces
  *every* fullerene to have exactly 12 pentagons, and specifically 20
  hexagons at n=60, regardless of how those faces are arranged, so a
  ring-count-only check wouldn't actually distinguish the one Ih-
  symmetric buckminsterfullerene isomer this module names from any of
  the many other, less symmetric C60 fullerene cages the same ring
  statistics also describe. A canonical-SMILES compare captures the
  exact connectivity instead.
- Deliberately not built on any of the existing ring machinery
  (`_aromatic.py`, `_polycyclic.py`, ...): none of those modules' scopes
  (ortho-fused chains, von Baeyer bridgeheads of degree <=3, etc.) come
  anywhere near a 12-pentagon/20-hexagon cage, and P-27's own general
  numbering system (the spiral algorithm, for substituted derivatives and
  other cage sizes) is a wholly different, novel piece of work out of
  scope here -- this module only recognizes the one retained name.

Explicitly out of scope: any substituent, any other fullerene cage size
or isomer (C70 etc.), and anything not exactly matching this one
structure. `has_fullerene_name` returns False for all of these, so
`core.py`'s existing dispatch continues to raise `UnsupportedStructure`
for them, unchanged.
"""

from rdkit import Chem

_FULLERENE_C60_SMILES = (
    "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1"
    "C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4="
    "C3C3=C5C(=C81)C%10=C23"
)
_FULLERENE_C60_CANONICAL = Chem.CanonSmiles(_FULLERENE_C60_SMILES)


def has_fullerene_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _FULLERENE_C60_CANONICAL


def name_fullerene(mol) -> str:
    return "[60]fullerene"
