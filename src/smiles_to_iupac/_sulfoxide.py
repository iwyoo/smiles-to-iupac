"""'sulfoxide' (P-63.6), implemented in `_chalcogen_oxide.py`."""

from ._chalcogen_oxide import ChalcogenOxide

_SULFOXIDE = ChalcogenOxide(
    atomic_num=16,
    element="sulfur",
    word="sulfoxide",
    prefix="sulfinyl",
    oxygens=1,
)

has_sulfoxide_shape = _SULFOXIDE.has_shape
name_sulfoxide = _SULFOXIDE.name
