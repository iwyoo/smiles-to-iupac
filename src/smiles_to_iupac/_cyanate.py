"""'R cyanate' (R-O-C#N), implemented in `_chalcocyanate.py`."""

from ._chalcocyanate import Cyanate

_CYANATE = Cyanate(
    atomic_num=8,
    symbol="O",
    element="oxygen",
    word="cyanate",
)

has_cyanate_shape = _CYANATE.has_shape
name_cyanate = _CYANATE.name
