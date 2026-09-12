"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
pentaphene, at the periphery bonds that don't touch a ring-fusion carbon
(P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf), via
`_fusion_locant_letter.py`'s shared matching engine.

`_LETTER_BY_PAIR` is derived by `_fusion_numbering.py`'s automatic
anchor-derivation path rather than hand-computed (PubChem CID 519935 as
reference): pentaphene's own mirror symmetry means no anchor set can
narrow its numbering to a single candidate (see that module's docstring),
so two known-compound anchors -- benzo[a]pentaphene (CID 21087903, letter
'a') and benzo[c]pentaphene (CID 123040, letter 'c') -- narrow it to
exactly the two mirror-image candidates; picking either names every
compound identically, so the tie-break below (lower reference-SMILES atom
index wins the 'a' bond) is an arbitrary but fixed internal convention,
verified to reproduce the two named anchors.

Of pentaphene's seven periphery bonds not touching a ring-fusion atom,
the resulting fused product at each was identified by building the
structure and matching InChIKey against PubChem:
- 'a'/'o' -- benzo[a]pentaphene (CID 21087903).
- 'b'/'n' -- hexaphene (CID 123042), confirmed by direct InChIKey match
  against PubChem's own "hexaphene" entry -- a retained name, cited as
  such rather than "benzo[b]pentaphene".
- 'c'/'m' -- benzo[c]pentaphene (CID 123040).
- 'h' -- fusing a ring here closes a bay region so that a pentaphene
  ring-fusion atom ends up shared by three rings (a branched/peri-fused
  arrangement, confirmed via direct construction), out of scope like
  every sibling module's peri-fused exclusion.

Scope, deliberately narrow, matching the sibling fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto pentaphene at the
'a', 'b', or 'c' bond. Explicitly out of scope (raise
`UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom, or the 'h' branched
  shape above.
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions
from ._fusion_numbering import letter_by_pair_candidates

_PENTAPHENE_REF = Chem.MolFromSmiles("C1=CC=C2C=C3C(=CC2=C1)C=CC4=CC5=CC=CC=C5C=C43")
# PubChem CID 21087903 (benzo[a]pentaphene, letter 'a') and CID 123040
# (benzo[c]pentaphene, letter 'c') as known-compound anchors.
_BENZO_A_PENTAPHENE = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=CC4=C(C=C32)C5=CC6=CC=CC=C6C=C5C=C4")
_BENZO_C_PENTAPHENE = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=C4C=CC5=CC6=CC=CC=C6C=C5C4=C3")
_NUMBERING_CANDIDATES = letter_by_pair_candidates(
    _PENTAPHENE_REF,
    [],
    compound_anchors=[(_BENZO_A_PENTAPHENE, "a"), (_BENZO_C_PENTAPHENE, "c")],
)
_LETTER_BY_PAIR = min(
    _NUMBERING_CANDIDATES,
    key=lambda letters: min(atom for pair, letter in letters.items() if letter == "a" for atom in pair),
)
_EXCLUDED_LETTERS = {"h"}
_RETAINED_NAME_BY_LETTER = {"b": "hexaphene"}

has_pentaphene_fusion_name, name_pentaphene_fusion = make_fusion_component_functions(
    "pentaphene",
    _PENTAPHENE_REF,
    _LETTER_BY_PAIR,
    _EXCLUDED_LETTERS,
    num_atoms=26,
    num_rings=6,
    retained_name_by_letter=_RETAINED_NAME_BY_LETTER,
)
