"""'ditelluride' and longer chalcogen chains R-Te(n)-R', implemented in `_dichalcogenide.py`."""

from ._dichalcogenide import Dichalcogenide

_DITELLURIDE = Dichalcogenide(
    atomic_num=52,
    symbol="Te",
    element="tellurium",
    word="ditelluride",
    base="telluride",
    stem="tellanyl",
    perol="pertellurol",
)

has_ditelluride_shape = _DITELLURIDE.has_shape
name_ditelluride = _DITELLURIDE.name
