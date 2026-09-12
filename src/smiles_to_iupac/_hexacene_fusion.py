"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
hexacene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
same mechanism as `_anthracene_fusion.py`/`_tetracene_fusion.py`/
`_pentacene_fusion.py`/`_phenanthrene_fusion.py`/`_pyrene_fusion.py`/
`_chrysene_fusion.py`/`_triphenylene_fusion.py`/`_fluoranthene_fusion.py`/
`_aceanthrylene_fusion.py`/`_acephenanthrylene_fusion.py` (all share
`_fusion_locant_letter.py`'s matching engine), this time for hexacene --
recognized directly by `_aromatic.py`'s general straight/angular ortho-
fused-chain algorithm (P-25.3.3.1.1) as the 6-ring linear case, the same
way tetracene/pentacene are -- as the base component. The eleventh base
overall, and the fourth (after anthracene/tetracene/pentacene) built on
a plain catacondensed chain.

Hexacene's own numbering (1,2,3,4,4a,5,5a,6,6a,7,7a,8,8a,9,10,11,12,12a,
13,13a,14,14a,15,15a,16,16a) was reconstructed from first principles
(SSSR + bond-ring-membership to find the periphery cycle, all 26 atoms
lying on it since hexacene, like the shorter acenes, is fully
catacondensed with no interior atom), then anchored against the only two
real methylhexacene structures PubChem indexes -- CID 149723163
(1-methyl) and CID 173022943 (2-methyl) -- via substructure match. Of
the six bonds not touching a ring-fusion atom (4a, 5a, 6a, 7a, 8a, 12a)
-- a, b, c, n, o, p -- hexacene's own symmetry (order 4, the same D2h
long-axis mirror + ring-swap rotation pair the shorter acenes have)
collapses these into just two distinct shapes: {a, c, n, p} (letter 'a'
wins per P-25.3.1.3's tie-break) and {b, o}. The {b, o} shape is the
linear extension -- the exact same compound (InChI-confirmed) as the
retained name 'heptacene' `_aromatic.py` recognizes directly (PR #552)
-- excluded here for the same reason `_anthracene_fusion.py`/
`_tetracene_fusion.py`/`_pentacene_fusion.py` exclude their own linear-
extension letters, so it falls through to that recognition unchanged.
The other, letter 'a', is a distinct, real registered compound:
'benzo[a]hexacene' (PubChem CID 23505249).

Scope, deliberately narrow, matching the other fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto hexacene at this
one position. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism), or the linear ('heptacene')
  shape.
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_HEXACENE_REF = Chem.MolFromSmiles("c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1")
_LETTER_BY_PAIR = {
    frozenset({25, 0}): "a",
    frozenset({0, 1}): "b",
    frozenset({1, 2}): "c",
    frozenset({12, 13}): "n",
    frozenset({13, 14}): "o",
    frozenset({14, 15}): "p",
}
_EXCLUDED_LETTERS = frozenset({"b"})

has_hexacene_fusion_name, name_hexacene_fusion = make_fusion_component_functions(
    "hexacene", _HEXACENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=30, num_rings=7
)
