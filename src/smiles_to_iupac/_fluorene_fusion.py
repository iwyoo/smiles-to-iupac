"""P-25.3.1.3 fusion-locant-letter citation for a plain benzo ring
ortho-fused onto fluorene, using `_fusion_locant_letter.py`'s engine.
Unlike `_fluoranthene_fusion.py`'s fully-aromatic base, fluorene's sp3
CH2 gets a fresh indicated-H locant on the fused system, per letter
hand-anchored against real PubChem CID 9195/9201/9150."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_locant_letter import find_fusion_letter

_FLUORENE_REF = Chem.MolFromSmiles("C1c2ccccc2-c2ccccc21")
_LETTER_BY_PAIR = {
    frozenset({2, 3}): "a",
    frozenset({3, 4}): "b",
    frozenset({4, 5}): "c",
    frozenset({8, 9}): "c",
    frozenset({9, 10}): "b",
    frozenset({10, 11}): "a",
}
_INDICATED_H_BY_LETTER = {"a": "11H", "b": "11H", "c": "7H"}


def _find_core(mol):
    return find_fusion_letter(
        mol,
        _FLUORENE_REF,
        _LETTER_BY_PAIR,
        frozenset(),
        num_atoms=17,
        num_rings=4,
        nonaromatic_ref_atoms=frozenset({0}),
    )


def has_fluorene_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_fluorene_fusion(mol) -> str:
    letter = _find_core(mol)
    if letter is None:
        raise UnsupportedStructure("this tetracyclic system is not a supported benzo-fused fluorene shape (see P-25.3.1.3)")
    return f"{_INDICATED_H_BY_LETTER[letter]}-benzo[{letter}]fluorene"
