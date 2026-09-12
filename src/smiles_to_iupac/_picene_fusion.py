"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
picene, at the periphery bonds that don't touch a ring-fusion carbon
(P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) -- picene
itself is already produced as a retained name by `_chrysene_fusion.py`
(chrysene + one ring, letter 'a'); this module is the next base further,
naming picene + one more ring (hexacyclic), via
`_fusion_locant_letter.py`'s shared matching engine.

`_LETTER_BY_PAIR` is derived by `_fusion_numbering.py`'s automatic
anchor-derivation path rather than hand-computed (PubChem CID 9162 as
reference): picene's own mirror symmetry means no anchor set can narrow
its numbering to a single candidate (see that module's docstring), so two
known-compound anchors -- benzo[b]picene (CID 123038, letter 'b') and
benzo[c]picene (CID 9168, letter 'c') -- narrow it to exactly the two
mirror-image candidates; picking either names every compound identically
(a mirror relabels every letter consistently), so the tie-break below
(lower reference-SMILES atom index wins the 'a' bond) is an arbitrary but
fixed internal convention, verified to reproduce the two named anchors.

Of picene's nine periphery bonds not touching a ring-fusion atom, the
resulting fused product at each was identified by building the structure
and matching InChIKey against PubChem:
- 'b'/'n' -- benzo[b]picene (CID 123038).
- 'c'/'m' -- benzo[c]picene (CID 9168).
- 'a'/'o' -- a real, registered compound (CID 13184858) with no trivial
  name indexed on PubChem; left out of scope here (adding systematic
  naming for it later doesn't touch this module's already-verified
  letters).
- 'f'/'j' and 's' -- fusing a ring here closes a bay region so that a
  picene ring-fusion atom ends up shared by three rings (a branched/
  peri-fused arrangement, confirmed via direct construction), a
  different citation mechanism, out of scope like every sibling module's
  peri-fused exclusion.

Scope, deliberately narrow, matching the sibling fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto picene at the 'b' or
'c' bond. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom, the 'a' shape above (no
  trivial name yet), or the 'f'/'s' branched shapes.
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions
from ._fusion_numbering import letter_by_pair_candidates

_PICENE_REF = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC=CC=C54")
# PubChem CID 123038 (benzo[b]picene, letter 'b') and CID 9168
# (benzo[c]picene, letter 'c') as known-compound anchors.
_BENZO_B_PICENE = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC6=CC=CC=C6C=C54")
_BENZO_C_PICENE = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=C4C=CC6=CC=CC=C65")
_NUMBERING_CANDIDATES = letter_by_pair_candidates(
    _PICENE_REF,
    [],
    compound_anchors=[(_BENZO_B_PICENE, "b"), (_BENZO_C_PICENE, "c")],
)
_LETTER_BY_PAIR = min(
    _NUMBERING_CANDIDATES,
    key=lambda letters: min(atom for pair, letter in letters.items() if letter == "a" for atom in pair),
)
_EXCLUDED_LETTERS = {"a", "f", "s"}

has_picene_fusion_name, name_picene_fusion = make_fusion_component_functions(
    "picene",
    _PICENE_REF,
    _LETTER_BY_PAIR,
    _EXCLUDED_LETTERS,
    num_atoms=26,
    num_rings=6,
)
