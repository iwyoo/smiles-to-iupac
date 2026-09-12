"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
heptacene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
same mechanism as `_anthracene_fusion.py`/`_tetracene_fusion.py`/
`_pentacene_fusion.py`/`_hexacene_fusion.py`/`_phenanthrene_fusion.py`/
`_pyrene_fusion.py`/`_chrysene_fusion.py`/`_triphenylene_fusion.py`/
`_fluoranthene_fusion.py`/`_aceanthrylene_fusion.py`/
`_acephenanthrylene_fusion.py` (all share `_fusion_locant_letter.py`'s
matching engine), this time for heptacene -- recognized directly by
`_aromatic.py`'s general straight/angular ortho-fused-chain algorithm
(P-25.3.3.1.1) as the 7-ring linear case, the same way tetracene/
pentacene/hexacene are -- as the base component. The twelfth base
overall, and the fifth (after anthracene/tetracene/pentacene/hexacene)
built on a plain catacondensed chain.

Heptacene's own numbering (1,2,3,4,4a,5,5a,6,6a,7,7a,8,8a,9,9a,10,11,12,
13,13a,14,14a,15,15a,16,16a,17,17a,18,18a) was reconstructed from first
principles (SSSR + bond-ring-membership to find the periphery cycle, all
34 atoms lying on it since heptacene, like the shorter acenes, is fully
catacondensed with no interior atom). Unlike the shorter acenes, PubChem
indexes no methylheptacene isomer at all to anchor the numbering
direction against -- but this doesn't leave the direction genuinely
ambiguous the way it would for an asymmetric base (see
`_aceanthrylene_fusion.py`'s docstring): heptacene's own automorphism
group (order 4, the same D2h long-axis mirror + ring-swap rotation pair
every shorter acene in this family has) already guarantees that whichever
of the 4 symmetric copies of "position 1" is chosen, the fusion-locant-
letter engine's own alphabetically-lowest tie-break picks the identical
winning letter for the identical physical shape -- so getting the
absolute numbering "backwards" relative to the (unpublished) official
numbering cannot change which name this module produces, only the
graph topology (periphery cycle, fusion-atom identification) needs to be
right, and that part is pure graph theory, independent of any external
anchor. Of the six bonds not touching a ring-fusion atom (4a, 5a, 6a,
7a, 8a, 9a, 13a) -- a, b, c, p, q, r (p/q/r being the middle three sides,
around positions 10-13) -- this symmetry collapses these into just two
distinct shapes: {a, c, p, r} (letter 'a' wins per P-25.3.1.3's tie-
break) and {b, q}. The {b, q} shape is the linear extension -- the exact
same compound (InChI-confirmed) as the retained name 'octacene'
`_aromatic.py` recognizes directly (PR #554) -- excluded here for the
same reason the shorter acenes' own modules exclude their linear-
extension letters, so it falls through to that recognition unchanged.
The other, letter 'a', is a real registered PubChem structure (CID
21087897, formula cross-checked as C34H20) but carries no name synonym
there to independently confirm 'benzo[a]heptacene' the way the shorter
acenes' equivalents did -- the name here rests on the fusion-locant-
letter algorithm itself (already validated against real named compounds
six times over in the sibling modules above) rather than an external
name match.

Scope, deliberately narrow, matching the other fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto heptacene at this
one position. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism), or the linear ('octacene')
  shape.
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_HEPTACENE_REF = Chem.MolFromSmiles("c1ccc2cc3cc4cc5cc6cc7ccccc7cc6cc5cc4cc3cc2c1")
_LETTER_BY_PAIR = {
    frozenset({29, 0}): "a",
    frozenset({0, 1}): "b",
    frozenset({1, 2}): "c",
    frozenset({14, 15}): "p",
    frozenset({15, 16}): "q",
    frozenset({16, 17}): "r",
}
_EXCLUDED_LETTERS = frozenset({"b"})

has_heptacene_fusion_name, name_heptacene_fusion = make_fusion_component_functions(
    "heptacene", _HEPTACENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=34, num_rings=8
)
