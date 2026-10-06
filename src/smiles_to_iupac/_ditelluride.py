"""'ditelluride' (R-Te-Te-R'), implemented in `_dichalcogenide.py`."""

from ._dichalcogenide import Dichalcogenide

_DITELLURIDE = Dichalcogenide(
    atomic_num=52,
    symbol="Te",
    element="tellurium",
    word="ditelluride",
    prefix="ditellanyl",
)

has_ditelluride_shape = _DITELLURIDE.has_shape
name_ditelluride = _DITELLURIDE.name
