"""Shared factory for the boilerplate every `_X_fusion.py` base-component
module repeats around `_fusion_locant_letter.find_fusion_letter`: a
`_find_core` wrapper plus a `has_..._fusion_name`/`name_..._fusion` pair
that raises `UnsupportedStructure` with a ring-count word derived from
`num_rings` (2026-09-12 refactor, no behavior change - see each caller's
own module docstring for the base-specific chemistry/citations/PubChem
CIDs, which stay in the calling module).

Each base module keeps its own chemistry data (reference SMILES,
`letter_by_pair`, `excluded_letters`, `num_atoms`, `num_rings`, optional
`nonaromatic_ref_atoms`/`retained_name_by_letter`) - only the wrapper
functions that consume that data move here.
"""

from ._common import UnsupportedStructure
from ._fusion_locant_letter import find_fusion_letter

_RING_COUNT_WORD = {
    3: "tricyclic",
    4: "tetracyclic",
    5: "pentacyclic",
    6: "hexacyclic",
    7: "heptacyclic",
    8: "octacyclic",
}


def make_fusion_component_functions(
    base_name,
    ref,
    letter_by_pair,
    excluded_letters,
    num_atoms,
    num_rings,
    nonaromatic_ref_atoms=frozenset(),
    retained_name_by_letter=None,
):
    """Returns `(has_fusion_name, name_fusion)` for one all-carbon
    retained-name aromatic base. `base_name`: used in the generated
    `benzo[letter]<base_name>` name and the `UnsupportedStructure`
    message. `retained_name_by_letter`: {letter -> retained name},
    for letters where the fused system already has its own retained
    name instead of the systematic `benzo[letter]...` form. The other
    parameters match `find_fusion_letter`'s own."""
    retained_name_by_letter = retained_name_by_letter or {}
    ring_count_word = _RING_COUNT_WORD[num_rings]

    def _find_core(mol):
        return find_fusion_letter(
            mol,
            ref,
            letter_by_pair,
            excluded_letters,
            num_atoms=num_atoms,
            num_rings=num_rings,
            nonaromatic_ref_atoms=nonaromatic_ref_atoms,
        )

    def has_fusion_name(mol) -> bool:
        return _find_core(mol) is not None

    def name_fusion(mol) -> str:
        letter = _find_core(mol)
        if letter is None:
            raise UnsupportedStructure(
                f"this {ring_count_word} system is not a supported benzo-fused "
                f"{base_name} shape (see P-25.3.1.3)"
            )
        retained_name = retained_name_by_letter.get(letter)
        if retained_name is not None:
            return retained_name
        return f"benzo[{letter}]{base_name}"

    return has_fusion_name, name_fusion
