"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
triphenylene, at the periphery bonds that don't touch a ring-fusion
carbon (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf), via
`_fusion_locant_letter.py`'s shared matching engine -- triphenylene
itself is `_branched_fused_aromatic.py`'s own retained name.

Triphenylene's numbering was pinned down the same empirical way as
chrysene's (`_chrysene_fusion.py`): 1-/2-methyltriphenylene (PubChem
CIDs 76131/618051) fixed the two distinct CH environments, and
triphenylene's D3h symmetry (order-6 automorphism group, confirmed via
`GetSubstructMatches`) collapses the nine periphery bonds not touching a
ring-fusion atom into just two distinct shapes: {a, c, g, i, m, o}
(min letter 'a') and {b, h, n} (min letter 'b').

Unlike every other base so far, triphenylene's three outer rings are
each fused to the *same* central ring rather than to each other in a
chain, so its "periphery" also passes through three bay-region bonds
(each directly connecting two ring-fusion atoms, e.g. atom pair (3, 17)
in the reference SMILES below) that belong to only one SSSR ring (the
central one) -- these still count as periphery per P-25.3.1.3 (each
faces the molecule's outer boundary, not another ring), but every
letter touching one of their endpoints is excluded the same way a bond
touching any other ring-fusion atom is.

Letter 'a's shape turns out to be the exact same compound as
`_chrysene_fusion.py`'s letter 'g' (confirmed by building the fused
structure and matching InChI), so only one of the two base components
may claim it -- P-25.3.2.4's parent-component seniority criteria decide
which. Chrysene and triphenylene tie on every criterion through (f) (same
ring count, same all-6-membered profile, no heteroatoms) and also on (g),
"the greatest number of rings in a horizontal row when drawn in the
preferred orientation" (P-25.3.2.3): `_fusion_orientation.py`'s algorithm
confirms both give 2. Criteria (h)/(i) (heteroatom locants) don't apply
either (no heteroatoms), so the real tiebreak is (j), "the lower locants
for the peripheral fusion carbon atoms" -- a numbering-based rule
(P-25.3.3.1) this project hasn't implemented yet. 'a' is excluded here so
chrysene claims the shape instead, which PubChem's own autonaming agrees
with independently, but is not yet a verified application of (j) itself;
revisit once a full peripheral-numbering algorithm exists. (An earlier
version of this comment guessed the tie extended to P-25.3.2.3.3's
quadrant sub-criteria (b)/(c)/(d) as well and fell through to plain
alphabetical order -- both halves of that guess were wrong: those
sub-criteria only pick a single preferred orientation for (g) itself, and
aren't separate P-25.3.2.4 list entries at all, and the actual next entry
after (g) is (j)'s locant rule, not alphabetical order. See
`_fusion_orientation.py`'s module docstring and issue #570.) Letter 'b'
doesn't match any retained name or other base's output; PubChem's own
autoname agrees independently ('benzo[b]triphenylene', CID 9164), so it's
supported here as the one genuinely new name this base adds.

Scope, deliberately narrow, matching the sibling fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto triphenylene at the
'b' bond. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (including a bay-region
  bond), or the 'a' shape above (already produced via chrysene).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_TRIPHENYLENE_REF = Chem.MolFromSmiles("c1ccc2c(c1)c1ccccc1c1ccccc21")
_LETTER_BY_PAIR = {
    frozenset({5, 0}): "a",
    frozenset({0, 1}): "b",
    frozenset({1, 2}): "c",
    frozenset({16, 15}): "g",
    frozenset({15, 14}): "h",
    frozenset({14, 13}): "i",
    frozenset({10, 9}): "m",
    frozenset({9, 8}): "n",
    frozenset({8, 7}): "o",
}
_EXCLUDED_LETTERS = {"a"}

has_triphenylene_fusion_name, name_triphenylene_fusion = make_fusion_component_functions(
    "triphenylene", _TRIPHENYLENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=22, num_rings=5
)
