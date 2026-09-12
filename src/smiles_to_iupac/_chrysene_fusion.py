"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
chrysene, at the periphery bonds that don't touch a ring-fusion carbon
(P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) -- chrysene
itself now supported as a retained name (`_phenanthrene_fusion.py`'s
letter 'a'), continuing the same pattern one base further, via
`_fusion_locant_letter.py`'s shared matching engine.

Chrysene's own standard numbering (1,2,3,4,4a,4b,5,6,6a,7,8,9,10,10a,10b,
11,12,12a) was pinned down empirically rather than from a Blue Book
numbering diagram (not extractable as text from the cached P2 chapter,
same issue noted in `general-fusion-naming.md`): six known-position
methylchrysenes (PubChem CIDs 18779-18782, 19427, 15564, i.e.
1-/2-/3-/4-/5-/6-methylchrysene) were each substructure-matched against
chrysene's own reference SMILES to fix which reference atom index is
which numbered position, then the periphery was walked from there and
cross-checked against chrysene's own automorphisms (order 2, confirmed
via `GetSubstructMatches(chrysene, ...)`) to find the C2 symmetry
(1<->7, 2<->8, 3<->9, 4<->10, 4a<->10a, 4b<->10b, 5<->11, 6<->12,
6a<->12a).

That gives P-25.3.1.3's continuous peripheral lettering as: a(1,2),
b(2,3), c(3,4), d(4,4a), e(4a,4b), f(4b,5), g(5,6), h(6,6a), i(6a,7),
j(7,8), k(8,9), l(9,10), m(10,10a), n(10a,10b), o(10b,11), p(11,12),
q(12,12a), r(12a,1). Of the eight bonds not touching a ring-fusion atom
(4a, 4b, 6a, 10a, 10b, 12a) -- a, b, c, g, j, k, l, p -- the C2 symmetry
above collapses these to four distinct shapes: {a, j} (picene, PubChem
CID 9162 -- confirmed by building the fused structure and matching
InChI against PubChem's own -- a retained name, cited as such rather
than as "benzo[a]chrysene"), {b, k}, {c, l}, and {g, p}. The latter
three don't match any retained name or existing PubChem common name
(each resolves, by direct PubChem structure lookup, only to a von
Baeyer-style systematic name -- CIDs 9163/9135/9140 respectively --
confirming they're real, registered compounds without their own
trivial name), so all three are supported here as
'benzo[b]chrysene'/'benzo[c]chrysene'/'benzo[g]chrysene'.

Scope, deliberately narrow, matching the sibling fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto chrysene at the
'a', 'b', 'c', or 'g' bond. Explicitly out of scope (raise
`UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system,
  different citation mechanism).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_locant_letter import find_fusion_letter

_CHRYSENE_REF = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43")
_LETTER_BY_PAIR = {
    frozenset({5, 0}): "a",
    frozenset({0, 1}): "b",
    frozenset({1, 2}): "c",
    frozenset({10, 11}): "g",
    frozenset({13, 14}): "j",
    frozenset({14, 15}): "k",
    frozenset({15, 16}): "l",
    frozenset({7, 6}): "p",
}
_EXCLUDED_LETTERS = set()
_RETAINED_NAME_BY_LETTER = {"a": "picene"}


def _find_core(mol):
    return find_fusion_letter(mol, _CHRYSENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=22, num_rings=5)


def has_chrysene_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_chrysene_fusion(mol) -> str:
    letter = _find_core(mol)
    if letter is None:
        raise UnsupportedStructure(
            "this pentacyclic system is not a supported benzo-fused "
            "chrysene shape (see P-25.3.1.3)"
        )
    retained_name = _RETAINED_NAME_BY_LETTER.get(letter)
    if retained_name is not None:
        return retained_name
    return f"benzo[{letter}]chrysene"
