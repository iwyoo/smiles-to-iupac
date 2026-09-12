"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
phenanthrene, at the periphery bonds that don't touch a ring-fusion
carbon (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
mirroring `_anthracene_fusion.py`'s identical mechanism for anthracene,
the other tricyclic all-carbon base component -- both now share
`_fusion_locant_letter.py`'s matching engine.

Phenanthrene's own numbering (P-25.3.3, 1,2,3,4,4a,4b,5,6,7,8,8a,9,10,10a)
gives P-25.3.1.3's continuous peripheral lettering as: a(1,2), b(2,3),
c(3,4), d(4,4a), e(4a,4b), f(4b,5), g(5,6), h(6,7), i(7,8), j(8,8a),
k(8a,9), l(9,10), m(10,10a), n(10a,1). Of the seven bonds not touching a
ring-fusion atom (4a, 4b, 8a, 10a) -- a, b, c, g, h, i, l -- phenanthrene's
own mirror symmetry (ring A <-> ring C) makes g/h/i the same compound as
c/b/a respectively (already found automatically by trying every
substructure-match automorphism, same as `_anthracene_fusion.py`), so
only three *distinct* shapes are possible: 'a' (chrysene, PubChem CID
9171 -- confirmed by building the fused structure and matching InChI --
a retained name, cited as such rather than as "benzo[a]phenanthrene"),
'b' (benzo[a]anthracene, CID 5954 -- the exact same compound
`_anthracene_fusion.py` already names via anthracene as the base
component instead, so excluded here to avoid a second, redundant route to
it), and 'l' (triphenylene, CID 9170 -- already a retained name,
`_branched_fused_aromatic.py`, excluded here for the same reason). So
this module claims 'chrysene' and the one genuinely new systematic name,
'benzo[c]phenanthrene' (PubChem CID 9136).

Scope, deliberately narrow, matching `_anthracene_fusion.py`: exactly one
plain, unsubstituted benzo ring ortho-fused onto phenanthrene at the 'a'
or 'c' bond. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system,
  different citation mechanism), or the 'b'/'l' shapes above (both
  already produced via a different, more senior base component).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_PHENANTHRENE_REF = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
_LETTER_BY_PAIR = {
    frozenset({1, 2}): "a",
    frozenset({0, 1}): "b",
    frozenset({0, 13}): "c",
    frozenset({9, 10}): "g",
    frozenset({8, 9}): "h",
    frozenset({7, 8}): "i",
    frozenset({4, 5}): "l",
}
_EXCLUDED_LETTERS = {"b", "l"}
_RETAINED_NAME_BY_LETTER = {"a": "chrysene"}

has_phenanthrene_fusion_name, name_phenanthrene_fusion = make_fusion_component_functions(
    "phenanthrene",
    _PHENANTHRENE_REF,
    _LETTER_BY_PAIR,
    _EXCLUDED_LETTERS,
    num_atoms=18,
    num_rings=4,
    retained_name_by_letter=_RETAINED_NAME_BY_LETTER,
)
