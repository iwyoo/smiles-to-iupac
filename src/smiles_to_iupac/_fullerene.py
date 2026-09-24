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

C70 is recognized the same way, as a second exact-match entry:

- Reference structure taken from PubChem CID 16131935 (InChIKey
  `ATLMFJTZZPOKLC-UHFFFAOYSA-N`)'s own connectivity SMILES, cross-checked
  here via RDKit: 70 all-carbon atoms, every atom degree 3, ring
  perception of exactly 12 five-membered and 25 six-membered rings.
- Named `(C70-D5h(6))[5,6]fullerene` per P-27.2.3: unlike C60 (which has
  only one fullerene isomer at that size, so P-27.2.2's trivial
  `[60]fullerene` form is unambiguous), C70 has multiple named isomers
  (e.g. D5h(5) is a distinct one), so the systematic PIN form carrying
  the point group is used instead of a bare `[70]fullerene` trivial form.

Explicitly out of scope: any substituent, any other fullerene cage size
or isomer, and anything not exactly matching one of these two
structures. `has_fullerene_name` returns False for all of these, so
`core.py`'s existing dispatch continues to raise `UnsupportedStructure`
for them, unchanged.
"""

from rdkit import Chem

_FULLERENE_C60_SMILES = (
    "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1"
    "C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4="
    "C3C3=C5C(=C81)C%10=C23"
)
_FULLERENE_C70_SMILES = (
    "C12=C3C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%10=C%10C8=C5C1=C%10C1=C%13C5=C8C1=C2"
    "C1=C3C2=C3C%10=C%13C%14=C3C1=C8C1=C3C5=C%12C5=C8C%11=C%11C9=C7C7=C9C6=C4C2=C2"
    "C%10=C4C(=C29)C2=C6C(=C8C8=C9C6=C4C%13=C9C(=C%141)C3=C85)C%11=C27"
)

_FULLERENE_NAMES = {
    Chem.CanonSmiles(_FULLERENE_C60_SMILES): "[60]fullerene",
    Chem.CanonSmiles(_FULLERENE_C70_SMILES): "(C70-D5h(6))[5,6]fullerene",
}


def has_fullerene_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _FULLERENE_NAMES


def name_fullerene(mol) -> str:
    return _FULLERENE_NAMES[Chem.MolToSmiles(mol)]
