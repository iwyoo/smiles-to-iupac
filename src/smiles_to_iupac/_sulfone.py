"""'sulfone' (P-63.6), implemented in `_chalcogen_oxide.py`."""

from ._chalcogen_oxide import ChalcogenOxide

_SULFONE = ChalcogenOxide(
    atomic_num=16,
    element="sulfur",
    word="sulfone",
    prefix="sulfonyl",
    oxygens=2,
)

has_sulfone_shape = _SULFONE.has_shape
name_sulfone = _SULFONE.name
