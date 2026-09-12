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
structure and matching InChI) -- chrysene is the senior base component
for that shape (already produced there), so 'a' is excluded here to
avoid a second, redundant route to it. Letter 'b' doesn't match any
retained name or other base's output; PubChem's own autoname agrees
independently ('benzo[b]triphenylene', CID 9164), so it's supported
here as the one genuinely new name this base adds.

Scope, deliberately narrow, matching the sibling fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto triphenylene at the
'b' bond. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (including a bay-region
  bond), or the 'a' shape above (already produced via chrysene).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_locant_letter import find_fusion_letter

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


def _find_core(mol):
    return find_fusion_letter(
        mol, _TRIPHENYLENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=22, num_rings=5
    )


def has_triphenylene_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_triphenylene_fusion(mol) -> str:
    letter = _find_core(mol)
    if letter is None:
        raise UnsupportedStructure(
            "this pentacyclic system is not a supported benzo-fused "
            "triphenylene shape (see P-25.3.1.3)"
        )
    return f"benzo[{letter}]triphenylene"
