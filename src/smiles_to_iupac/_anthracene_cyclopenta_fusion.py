"""P-25.3.1.3 fusion-locant-letter citation for a cyclopenta ring
ortho-fused onto anthracene, via `_fusion_locant_letter.py`'s engine
generalized to a non-aromatic 3-atom attachment, plus
`_fusion_numbering_general.py`'s fresh whole-system indicated-H locant.
Anchored against real PubChem CID 21087838/57983840/19845433."""

from rdkit import Chem

from ._anthracene_fusion import _ANTHRACENE_REF, _LETTER_BY_PAIR
from ._common import UnsupportedStructure
from ._fusion_locant_letter import find_fusion_letter
from ._fusion_numbering_general import general_peripheral_numbering

_NUM_ATOMS = 17
_NUM_RINGS = 4


def _find_letter(mol):
    return find_fusion_letter(
        mol,
        _ANTHRACENE_REF,
        _LETTER_BY_PAIR,
        frozenset(),
        num_atoms=_NUM_ATOMS,
        num_rings=_NUM_RINGS,
        num_extra_atoms=3,
        nonaromatic_extra_atom_count=3,
    )


def has_anthracene_cyclopenta_fusion_name(mol) -> bool:
    return _find_letter(mol) is not None


def name_anthracene_cyclopenta_fusion(mol) -> str:
    letter = _find_letter(mol)
    if letter is None:
        raise UnsupportedStructure(
            "this tetracyclic system is not a supported cyclopenta-fused "
            "anthracene shape (see P-25.3.1.3)"
        )
    numbering = general_peripheral_numbering(mol)
    indicated_h = next(
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 2
    )
    return f"{numbering[indicated_h]}H-cyclopenta[{letter}]anthracene"
