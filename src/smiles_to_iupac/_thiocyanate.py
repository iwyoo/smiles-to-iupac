"""'R thiocyanate' (R-S-C#N), implemented in `_chalcocyanate.py`."""

from ._chalcocyanate import Cyanate

_THIOCYANATE = Cyanate(
    atomic_num=16,
    symbol="S",
    element="sulfur",
    word="thiocyanate",
)

has_thiocyanate_shape = _THIOCYANATE.has_shape
name_thiocyanate = _THIOCYANATE.name
