"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
pentacene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
same mechanism as `_anthracene_fusion.py`/`_tetracene_fusion.py`/
`_phenanthrene_fusion.py`/`_pyrene_fusion.py`/`_chrysene_fusion.py`/
`_triphenylene_fusion.py`/`_fluoranthene_fusion.py`/
`_aceanthrylene_fusion.py`/`_acephenanthrylene_fusion.py` (all share
`_fusion_locant_letter.py`'s matching engine), this time for pentacene --
already recognized directly by `_aromatic.py`'s general straight/angular
ortho-fused-chain algorithm (P-25.3.3.1.1) as the 5-ring linear case, the
same way tetracene is -- as the base component. The tenth base overall,
and the third (after anthracene/tetracene) built on a plain
catacondensed chain.

Pentacene's own numbering (1,2,3,4,4a,5,5a,6,6a,7,7a,8,9,10,11,11a,12,
12a,13,13a,14,14a) was reconstructed from first principles (SSSR + bond-
ring-membership to find the periphery cycle, all 22 atoms lying on it
since pentacene, like anthracene/tetracene, is fully catacondensed with
no interior atom), then anchored against the only two real
methylpentacene structures PubChem indexes -- CID 23505239 (1-methyl)
and CID 5246912 (2,3,9,10-tetramethyl) -- via substructure match. This
gives P-25.3.1.3's continuous peripheral lettering as: a(1,2), b(2,3),
c(3,4), d(4,4a), e(4a,5), f(5,5a), g(5a,6), h(6,6a), i(6a,7), j(7,7a),
k(7a,8), l(8,9), m(9,10), n(10,11), o(11,11a), p(11a,12), q(12,12a),
r(12a,13), s(13,13a), t(13a,14), u(14,14a), v(14a,1). Of the six bonds
not touching a ring-fusion atom (4a, 5a, 6a, 7a, 11a, 12a, 13a, 14a) --
a, b, c, l, m, n -- pentacene's own symmetry (order 4, the same D2h
long-axis mirror + ring-swap rotation pair anthracene/tetracene have)
collapses these into just two distinct shapes: {a, c, l, n} (letter 'a'
wins per P-25.3.1.3's tie-break) and {b, m}. The {b, m} shape is the
linear extension -- the exact same compound (InChI-confirmed) as the
retained name 'hexacene', which `_aromatic.py`'s retained-name table now
covers directly -- excluded here for the same reason
`_anthracene_fusion.py`/`_tetracene_fusion.py` exclude their own
linear-extension letters, so it falls through to that recognition
unchanged. The other, letter 'a', is a distinct, real registered
compound: 'benzo[a]pentacene' (PubChem CID 67482, CAS 239-98-5, also
known as "isohexaphene").

Scope, deliberately narrow, matching the other fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto pentacene at this
one position. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism), or the linear ('hexacene')
  shape.
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_PENTACENE_REF = Chem.MolFromSmiles("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1")
_LETTER_BY_PAIR = {
    frozenset({21, 0}): "a",
    frozenset({0, 1}): "b",
    frozenset({1, 2}): "c",
    frozenset({10, 11}): "l",
    frozenset({11, 12}): "m",
    frozenset({12, 13}): "n",
}
_EXCLUDED_LETTERS = frozenset({"b"})

has_pentacene_fusion_name, name_pentacene_fusion = make_fusion_component_functions(
    "pentacene", _PENTACENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=26, num_rings=6
)
