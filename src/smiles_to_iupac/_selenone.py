"""'selenone' (P-63.6), implemented in `_chalcogen_oxide.py`."""

from ._chalcogen_oxide import ChalcogenOxide

_SELENONE = ChalcogenOxide(
    atomic_num=34,
    element="selenium",
    word="selenone",
    prefix="selenonyl",
    oxygens=2,
)

has_selenone_shape = _SELENONE.has_shape
name_selenone = _SELENONE.name
