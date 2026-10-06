"""'R-N=C=Se' (isoselenocyanate, P-61.8), implemented in `_iso_cyanate.py`."""

from ._iso_cyanate import Isocyanate

_ISOSELENOCYANATE = Isocyanate(
    atomic_num=34,
    symbol="Se",
    element="selenium",
    infix="seleno",
)

has_isoselenocyanate_shape = _ISOSELENOCYANATE.has_shape
name_isoselenocyanate = _ISOSELENOCYANATE.name
