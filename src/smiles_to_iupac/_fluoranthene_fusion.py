"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
fluoranthene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
same mechanism as `_anthracene_fusion.py`/`_phenanthrene_fusion.py`/
`_pyrene_fusion.py`/`_chrysene_fusion.py`/`_triphenylene_fusion.py` (all
share `_fusion_locant_letter.py`'s matching engine), this time for
fluoranthene (`_peri_fused_aromatic.py`'s own retained-name reference
SMILES) as the base component -- fluoranthene is itself peri-fused (one
interior atom, like pyrene's two), so this is the second peri-fused base.

Fluoranthene's own numbering (1,2,3,3a,4,5,6,6a,6b,7,8,9,10,10a,10b,10c,
with 10c the sole interior atom) was reconstructed here from first
principles (SSSR + bond-ring-membership to find the periphery cycle and
the one interior atom, then anchored against five real
monomethylfluoranthene PubChem structures -- CID 146998/35702/15565/
146751/146538 for the 1-/2-/3-/7-/8-methyl isomers respectively -- via
substructure match, which fixed both the numbering direction and which
of fluoranthene's two mirror-symmetric branches is which) gives
P-25.3.1.3's continuous peripheral lettering as: a(1,2), b(2,3), c(3,3a),
d(3a,4), e(4,5), f(5,6), g(6,6a), h(6a,6b), i(6b,7), j(7,8), k(8,9),
l(9,10), m(10,10a), n(10a,10b), o(10b,1). Of the seven bonds not touching
a ring-fusion atom (3a, 6a, 6b, 10a, 10b; the interior atom 10c never
touches the periphery at all) -- a, b, e, f, j, k, l -- fluoranthene's
own mirror symmetry (its only nontrivial automorphism, confirmed by
trying every substructure-match automorphism, same method as the other
fusion modules) collapses these seven into four distinct shapes:
{a, f} = 'benzo[a]fluoranthene' (PubChem CID 9146), {b, e} =
'benzo[b]fluoranthene' (CID 9153), {j, l} = 'benzo[j]fluoranthene' (CID
9152 -- PubChem resolves the name 'benzo[l]fluoranthene' to this same
CID), and {k} alone (self-paired by the symmetry) = 'benzo[k]fluoranthene'
(CID 9158). Each of the four was independently verified by building the
fused structure directly and comparing its InChI against the real
compound's PubChem InChI -- all four match exactly, and (unlike the
anthracene/phenanthrene/chrysene cases) none collides with an
already-recognized retained name, so all four are supported here.

Scope, deliberately narrow, matching the other fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto fluoranthene at one
of these four distinct positions. Explicitly out of scope (raise
`UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_FLUORANTHENE_REF = Chem.MolFromSmiles("C1=CC=C2C(=C1)C3=CC=CC4=C3C2=CC=C4")
_LETTER_BY_PAIR = {
    frozenset({13, 14}): "a",
    frozenset({14, 15}): "b",
    frozenset({8, 9}): "e",
    frozenset({7, 8}): "f",
    frozenset({0, 5}): "j",
    frozenset({0, 1}): "k",
    frozenset({1, 2}): "l",
}
_EXCLUDED_LETTERS = frozenset()

has_fluoranthene_fusion_name, name_fluoranthene_fusion = make_fusion_component_functions(
    "fluoranthene", _FLUORANTHENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=20, num_rings=5
)
