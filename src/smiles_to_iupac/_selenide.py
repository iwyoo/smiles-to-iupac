"""'selenide' (R-Se-R'), implemented in `_chalcogenide.py`."""

from ._chalcogenide import Chalcogenide

_SELENIDE = Chalcogenide(
    atomic_num=34,
    symbol="Se",
    element="selenium",
    word="selenide",
    prefix="selanyl",
)

has_selenide_shape = _SELENIDE.has_shape
name_selenide = _SELENIDE.name
