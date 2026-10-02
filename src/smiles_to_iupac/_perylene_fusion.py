"""P-25.3.1.3 fusion-locant-letter naming for a plain benzo ring ortho-
fused onto perylene, reusing perylene's own peripheral numbering
unchanged (the result stays fully mancude). Real structures exist for
letters 'a' (CID 115240) and 'b' (CID 67455), though PubChem's own PIN
software uses von Baeyer numbering -- verified here by rule derivation."""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_PERYLENE_REF = Chem.MolFromSmiles("C1=CC2=C3C(=C1)C4=CC=CC5=C4C(=CC=C5)C3=CC=C2")
_LETTER_BY_PAIR = {
    frozenset({13, 14}): "a",
    frozenset({14, 15}): "b",
    frozenset({15, 10}): "c",
    frozenset({10, 9}): "d",
    frozenset({9, 8}): "e",
    frozenset({8, 7}): "f",
    frozenset({7, 6}): "g",
    frozenset({6, 4}): "h",
    frozenset({4, 5}): "i",
    frozenset({5, 0}): "j",
    frozenset({0, 1}): "k",
    frozenset({1, 2}): "l",
    frozenset({2, 19}): "m",
    frozenset({19, 18}): "n",
    frozenset({18, 17}): "o",
    frozenset({17, 16}): "p",
    frozenset({16, 12}): "q",
    frozenset({12, 13}): "r",
}

has_perylene_fusion_name, name_perylene_fusion = make_fusion_component_functions(
    "perylene", _PERYLENE_REF, _LETTER_BY_PAIR, frozenset(), num_atoms=24, num_rings=6
)
