"""'R-N=C=O' (isocyanate, P-61.8), implemented in `_iso_cyanate.py`."""

from ._iso_cyanate import Isocyanate

_ISOCYANATE = Isocyanate(
    atomic_num=8,
    symbol="O",
    element="oxygen",
    infix="",
)

has_isocyanate_shape = _ISOCYANATE.has_shape
name_isocyanate = _ISOCYANATE.name
