"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
tetracene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
same mechanism as `_anthracene_fusion.py`/`_phenanthrene_fusion.py`/
`_pyrene_fusion.py`/`_chrysene_fusion.py`/`_triphenylene_fusion.py`/
`_fluoranthene_fusion.py`/`_aceanthrylene_fusion.py`/
`_acephenanthrylene_fusion.py` (all share `_fusion_locant_letter.py`'s
matching engine), this time for tetracene -- already recognized directly
by `_aromatic.py`'s general straight/angular ortho-fused-chain algorithm
(P-25.3.3.1.1) as the 4-ring linear case, the same way anthracene is --
as the base component. The ninth base overall, and the second (after
anthracene) built on a plain catacondensed chain rather than a
peri-fused or already-dedicated-module retained name.

Tetracene's own numbering (1,2,3,4,4a,5,5a,6,6a,7,8,9,10,10a,11,11a,12,
12a) was reconstructed from first principles (SSSR + bond-ring-
membership to find the periphery cycle, all 18 atoms lying on it since
tetracene, like anthracene, is fully catacondensed with no interior
atom), then anchored against real monomethyltetracene/monomethyl-
naphthacene PubChem structures -- CID 23505253/15660752/23505242 for the
1-/2-/5-methyl isomers respectively -- via substructure match. This
gives P-25.3.1.3's continuous peripheral lettering as: a(1,2), b(2,3),
c(3,4), d(4,4a), e(4a,5), f(5,5a), g(5a,6), h(6,6a), i(6a,7), j(7,8),
k(8,9), l(9,10), m(10,10a), n(10a,11), o(11,11a), p(11a,12), q(12,12a),
r(12a,1). Of the six bonds not touching a ring-fusion atom (4a, 5a, 6a,
10a, 11a, 12a) -- a, b, c, j, k, l -- tetracene's own symmetry (order 4,
the same D2h long-axis mirror + ring-swap rotation pair anthracene has)
collapses these into just two distinct shapes: {a, c, j, l} (letter 'a'
wins per P-25.3.1.3's tie-break) and {b, k}. The {b, k} shape is the
linear extension -- the exact same compound (InChI-confirmed) as the
plain retained name 'pentacene' `_aromatic.py` already recognizes
directly -- and is excluded here for the same reason
`_anthracene_fusion.py` excludes its own letter 'b' (the tetracene
shape). The other, letter 'a', is a distinct, real registered compound:
'benzo[a]tetracene' (PubChem CID 67470, and PubChem's own IUPACName
property agrees exactly).

Scope, deliberately narrow, matching the other fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto tetracene at this
one position. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism), or the linear ('pentacene')
  shape.
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_TETRACENE_REF = Chem.MolFromSmiles("c1ccc2cc3cc4ccccc4cc3cc2c1")
_LETTER_BY_PAIR = {
    frozenset({1, 2}): "a",
    frozenset({0, 1}): "b",
    frozenset({0, 17}): "c",
    frozenset({10, 11}): "j",
    frozenset({9, 10}): "k",
    frozenset({8, 9}): "l",
}
_EXCLUDED_LETTERS = frozenset({"b"})

has_tetracene_fusion_name, name_tetracene_fusion = make_fusion_component_functions(
    "tetracene", _TETRACENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=22, num_rings=5
)
