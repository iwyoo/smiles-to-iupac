"""'sulfide' (R-S-R'), implemented in `_chalcogenide.py`."""

from ._chalcogenide import Chalcogenide

_SULFIDE = Chalcogenide(
    atomic_num=16,
    symbol="S",
    element="sulfur",
    word="sulfide",
    prefix="sulfanyl",
)

has_sulfide_shape = _SULFIDE.has_shape
name_sulfide = _SULFIDE.name
