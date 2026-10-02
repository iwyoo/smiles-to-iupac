"""P-25.3.1.3 fusion-locant-letter citation for a plain benzo ring
ortho-fused onto azulene, using `_fusion_locant_letter.py`'s engine.
Azulene's mirror symmetry collapses six candidate bonds into three
classes; two are anchored against real PubChem CID 601271/14150596,
the third ({f,g}) is structurally identical by the same symmetry."""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_AZULENE_REF = Chem.MolFromSmiles("C1=CC2=CC=CC=CC2=C1")
_LETTER_BY_PAIR = {
    frozenset({9, 0}): "a",
    frozenset({0, 1}): "a",
    frozenset({3, 4}): "e",
    frozenset({6, 7}): "e",
    frozenset({4, 5}): "f",
    frozenset({5, 6}): "f",
}
_EXCLUDED_LETTERS = frozenset()

has_azulene_fusion_name, name_azulene_fusion = make_fusion_component_functions(
    "azulene",
    _AZULENE_REF,
    _LETTER_BY_PAIR,
    _EXCLUDED_LETTERS,
    num_atoms=14,
    num_rings=3,
)
