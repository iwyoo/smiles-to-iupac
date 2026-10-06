"""'selenoxide' (P-63.6), implemented in `_chalcogen_oxide.py`."""

from ._chalcogen_oxide import ChalcogenOxide

_SELENOXIDE = ChalcogenOxide(
    atomic_num=34,
    element="selenium",
    word="selenoxide",
    prefix="seleninyl",
    oxygens=1,
)

has_selenoxide_shape = _SELENOXIDE.has_shape
name_selenoxide = _SELENOXIDE.name
