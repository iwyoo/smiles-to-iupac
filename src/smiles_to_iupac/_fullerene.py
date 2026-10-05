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

C76 is recognized the same way, as a third exact-match entry:

- Reference structure taken from PubChem CID 56846604 (InChIKey
  `DEJYFPHYOINFQD-UHFFFAOYSA-N`)'s own connectivity SMILES, cross-checked
  here via RDKit: 76 all-carbon atoms, every atom degree 3, ring
  perception of exactly 12 five-membered and 28 six-membered rings.
- Named `(C76-D2)[5,6]fullerene`: D2 is C76's sole isolable
  isolated-pentagon-rule isomer, so it's the isomer PubChem/CAS register
  under the plain "Fullerene C76" name -- the same point-group-qualified
  PIN form as C70 above, since C76 (like C70) has more than one distinct
  fullerene isomer at that atom count.

A single skeletal-silicon replacement on the C60-Ih cage is recognized
the same way, as a fourth exact-match entry:

- Reference structure taken from PubChem CID 101063510 (formula C59Si)'s
  own connectivity SMILES, cross-checked here via RDKit: 60 skeletal
  atoms total, every carbon degree 3, the one silicon atom degree 3 with
  0 H and 0 radical electrons (silicon's standard valence 4 is satisfied
  the same way carbon's is in the cage's aromatic system, unlike a
  trivalent heteroatom which would need an added hydrogen), ring
  perception of exactly 12 five-membered and 20 six-membered rings --
  identical to the plain C60-Ih entry above, just one vertex relabeled.
- Named `sila(C60-Ih)[5,6]fullerene` per P-27.5.1's own worked example,
  with no locant for the replacement position: every one of C60-Ih's 60
  carbon atoms is symmetry-equivalent (a single orbit under the Ih point
  group), so no locant is needed to disambiguate regardless of which
  vertex was replaced -- the point-group symbol in the name is the
  parent fullerene's own symmetry label, unchanged by the substitution.
- A trivalent heteroatom replacement (aza, phospha, ...), a
  multi-heteroatom replacement, or any heteroatom replacement on the
  C70/C76 cages (both have multiple symmetry-inequivalent carbon orbits,
  so a single replacement there would need a locant, which needs real
  fullerene numbering (P-27.3) this project doesn't have yet) are each
  out of scope.

A single methylene fused across an intact 6,6-bond of the C60-Ih cage
(P-27.6.1's ortho-fused-cyclopropane case -- not P-27.4.1's
homofullerene, see below) is recognized the same way, as a fifth
exact-match entry:

- Reference structure taken from PubChem CID 11422743 (formula C61H2,
  InChIKey `JURXXEICUOUOOF-UHFFFAOYSA-N`)'s own connectivity SMILES,
  cross-checked here via RDKit: 61 skeletal atoms, ring perception of
  exactly 12 five-membered, 20 six-membered, and one new three-membered
  ring. The two former cage carbons at the base of that triangle are
  *still bonded to each other* (each has degree 4, one neighbor being
  the other) in addition to each bonding to the new CH2 -- a fused
  cyclopropane sharing that intact 6,6-bond as its fusion edge, not a
  bond broken and replaced by the methylene. The two bridgehead atoms
  both belong to two six-membered rings together (in addition to the
  new cyclopropane ring), confirming the fusion sits on a 6,6-bond
  (between two hexagons), not a 5,6-bond.
- Named `3'H-cyclopropa[1,9](C60-Ih)[5,6]fullerene` -- P-27.6.1's own
  literal worked example (`tmp/bluebook/P2.txt` ~7617, diagram on PDF
  p.147: atoms 1 and 9 are drawn directly bonded to each other, both
  also bonded to the new ring carbon 3'). This is a genuinely different
  molecule from P-27.4.1's homofullerene
  (`1(9)aH-1(9)a-homo(C60-Ih)[5,6]fullerene`, worked-example diagram on
  PDF p.143): there the methylene is inserted *into* the 1-9 bond, so
  1 and 9 end up bonded only through the new CH2 and not to each other
  -- a ring-expansion with no cyclopropane, not the structure this
  reference SMILES encodes (an earlier version of this module conflated
  the two and misnamed this PubChem structure as the homofullerene).
  C60-Ih's 6,6-bonds form a single symmetry orbit (all 30 of them
  equivalent under the Ih point group, the well-known basis for "the"
  6,6-bond in fullerene chemistry, e.g. PCBM-style methanofullerenes),
  so P-27.3's systematic numbering always normalizes a 6,6-bond fusion
  to locants 1/9 regardless of which physical bond was used -- the
  literal worked-example locant string applies to any such structure,
  not just one specific orientation, so no separate general numbering
  implementation is needed for this one exact-match case.
- A fusion across a 5,6-bond instead, more than one such fusion, the
  true P-27.4.1 homofullerene (no registered PubChem structure found
  for it), the same cyclopropane fusion on the C70/C76/sila-C60 cages,
  and the nor-/seco-/cyclo- prefix families (P-27.4.2/.3/.4) are each
  out of scope -- separate, unresearched follow-up work.

C84 cage sizes are handled differently from the five exact-match entries
above: C84 has 24 distinct isolated-pentagon-rule isomers (unlike C60/C70/
C76, which each have only one isomer -- or one isolable isomer -- at that
size), so no single reference SMILES can cover it. `_fullerene_spiral.py`
implements the real fix for what issue #997 called a missing "isomer-
atlas/spiral-code cross-reference tool": it computes the given cage's own
canonical ring spiral (the Fowler-Manolopoulos algorithm, from the
structure's planar embedding, not from any lookup shortcut) and matches it
against the 24 known IPR isomers' published spirals. See that module's
docstring for the full algorithm and its sourcing/validation.
`has_fullerene_name`/`name_fullerene` fall through to it for any
84-carbon cage that isn't one of this module's own five hardcoded
entries.

A substituted C60/C70 cage (as parent with hydro prefixes, or as a group with free valences and added hydrogen,
P-29.3.4.1, P-6) is named through `_fullerene_numbering.py`; no other substituted cage is.

Explicitly out of scope: any other substituent, any heteroatom replacement or
cyclopropane fusion other than the two single-site C60 cases above, any
cage size other than 60/70/76/84, and (even at size 84) any non-IPR
isomer or isomer not among the 24 in `_fullerene_spiral.py`'s table.
`has_fullerene_name` returns False for all of these, so `core.py`'s
existing dispatch continues to raise `UnsupportedStructure` for them,
unchanged.
"""

from functools import lru_cache

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure
from ._fullerene_spiral import match_c84_isomer

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
_FULLERENE_C76_SMILES = (
    "C12=C3C4=C5C6=C1C7=C8C2=C9C1=C2C%10=C%11C(=C13)C1=C4C3=C4C%12=C%13C3=C5C3=C6"
    "C5=C7C6=C7C5=C5C3=C%13C3=C%13C%12=C%12C%14=C%15C%16=C%17C%18=C%19C%20=C%16C"
    "(=C%13%14)C(=C35)C7=C%20C3=C%19C5=C7C%18=C%13C%17=C%14C%15=C%15C%12=C4C1=C%11"
    "C%15=C%14C%10=C%13C2=C7C9=C5C8=C63"
)

_SILA_C60_SMILES = (
    "C12=C3C4=C5C6=C7C8=C9C(=C61)C1=C6C%10=C%11C(=C21)C1=C3C2=C4C3=C4C5=C7C5=C7C8="
    "C8C9=C6C6=C%10C9=C%10C%11=C1C1=C2C2=C3C3=C%11C%12=C2C1=C%10C1=C9C2=C9C(=C7C(="
    "C%11[Si]9=C1%12)C5=C43)C8=C62"
)

_CYCLOPROPA_C60_SMILES = (
    "C1C23C14C5=C6C7=C8C9=C1C%10=C%11C%12=C%13C%14=C%10C%10=C1C1=C%15C%16=C%17C%18="
    "C%19C%20=C%21C%22=C%23C%24=C%25C%26=C(C7=C9C%11=C%26C%12=C%24C%22=C%13C%20=C%14"
    "C%18=C%10%16)C7=C%25C9=C(C4=C76)C4=C2C(=C%17C3=C%15C5=C81)C%19=C%21C4=C%239"
)

# Reference structure for one of C84's 24 IPR isomers (isomer #24, D6h),
# used only as a test fixture for `_fullerene_spiral.py`'s general spiral-
# matching mechanism -- see that module's docstring. Taken from PubChem CID
# 133108900 (InChIKey FQRWAZOLUJHNDT-UHFFFAOYSA-N).
_C84_D6H_ISOMER_24_SMILES = (
    "C12=C3C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%14=C%15C%16=C%17C%18=C%19C%20=C%21C"
    "%22=C(C1=C1C(=C36)C(=C7%10)C(=C%11%14)C(=C%18%15)C1=C%22%19)C1=C2C2=C4C3=C4"
    "C5=C8C5=C6C9=C%12C7=C8C%13=C%16C9=C%10C%17=C%20C%11=C%12C%21=C1C1=C2C2=C3C3="
    "C%13C%14=C2C1=C%12C1=C%14C2=C(C%10=C%111)C9=C8C1=C2C%13=C(C5=C43)C6=C71"
)

_FULLERENE_NAMES = {
    Chem.CanonSmiles(_FULLERENE_C60_SMILES): "[60]fullerene",
    Chem.CanonSmiles(_FULLERENE_C70_SMILES): "(C70-D5h(6))[5,6]fullerene",
    Chem.CanonSmiles(_FULLERENE_C76_SMILES): "(C76-D2)[5,6]fullerene",
    Chem.CanonSmiles(_SILA_C60_SMILES): "sila(C60-Ih)[5,6]fullerene",
    Chem.CanonSmiles(_CYCLOPROPA_C60_SMILES): "3'H-cyclopropa[1,9](C60-Ih)[5,6]fullerene",
}


def has_fullerene_name(mol) -> bool:
    if Chem.MolToSmiles(mol) in _FULLERENE_NAMES:
        return True
    return match_c84_isomer(mol) is not None


def name_fullerene(mol) -> str:
    smi = Chem.MolToSmiles(mol)
    if smi in _FULLERENE_NAMES:
        return _FULLERENE_NAMES[smi]
    return match_c84_isomer(mol)


def is_fullerene_cage(mol, atoms) -> bool:
    """`atoms` form a closed carbon cage of twelve five-membered rings and only five- or six-membered rings."""
    atoms = set(atoms)
    if len(atoms) < 20 or any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in atoms):
        return False
    rings = [r for r in mol.GetRingInfo().AtomRings() if set(r) <= atoms]
    return sum(len(r) == 5 for r in rings) == 12 and all(len(r) in (5, 6) for r in rings)


def has_substituted_fullerene_cage(mol) -> bool:
    ring_info = mol.GetRingInfo()
    if sum(len(r) == 5 for r in ring_info.AtomRings()) < 12:
        return False
    cage = {a for r in ring_info.AtomRings() for a in r}
    return mol.GetNumAtoms() > len(cage) and is_fullerene_cage(mol, cage)


def _cage_graph(mol, atoms):
    atoms = set(atoms)
    return nx.Graph(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
    )


@lru_cache(maxsize=None)
def _reference_graph(smiles):
    return _cage_graph(Chem.MolFromSmiles(smiles), range(Chem.MolFromSmiles(smiles).GetNumAtoms()))


_NUMBERED_CAGE_STEMS = (
    (_FULLERENE_C60_SMILES, "(C60-Ih)[5,6]fullerene"),
    (_FULLERENE_C70_SMILES, "(C70-D5h(6))[5,6]fullerene"),
)


def numbered_cage_stem(mol, atoms):
    """The substitutive stem of the C60-Ih or C70-D5h(6) cage formed by `atoms`, or None for any other cage."""
    graph = _cage_graph(mol, atoms)
    for smiles, stem in _NUMBERED_CAGE_STEMS:
        reference = _reference_graph(smiles)
        if graph.number_of_nodes() == reference.number_of_nodes() and nx.is_isomorphic(graph, reference):
            return stem
    return None


def require_defined_fullerene_numbering(mol, atoms):
    if is_fullerene_cage(mol, atoms) and numbered_cage_stem(mol, atoms) is None:
        raise UnsupportedStructure(
            "fullerene locants cannot be derived: Fu-3.1 states its rules suffice only for the C60-Ih and "
            "C70-D5h(6) cages and that more rules are needed for other fullerenes"
        )
