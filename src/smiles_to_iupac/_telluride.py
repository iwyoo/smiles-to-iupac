"""'telluride' (R-Te-R'), implemented in `_chalcogenide.py`."""

from ._chalcogenide import Chalcogenide

_TELLURIDE = Chalcogenide(
    atomic_num=52,
    symbol="Te",
    element="tellurium",
    word="telluride",
    prefix="tellanyl",
)

has_telluride_shape = _TELLURIDE.has_shape
name_telluride = _TELLURIDE.name
