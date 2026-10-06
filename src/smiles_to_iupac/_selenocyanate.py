"""'R selenocyanate' (R-Se-C#N), implemented in `_chalcocyanate.py`."""

from ._chalcocyanate import Cyanate

_SELENOCYANATE = Cyanate(
    atomic_num=34,
    symbol="Se",
    element="selenium",
    word="selenocyanate",
)

has_selenocyanate_shape = _SELENOCYANATE.has_shape
name_selenocyanate = _SELENOCYANATE.name
