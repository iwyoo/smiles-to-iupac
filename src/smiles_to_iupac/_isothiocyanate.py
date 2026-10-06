"""'R-N=C=S' (isothiocyanate, P-61.8), implemented in `_iso_cyanate.py`."""

from ._iso_cyanate import Isocyanate

_ISOTHIOCYANATE = Isocyanate(
    atomic_num=16,
    symbol="S",
    element="sulfur",
    infix="thio",
)

has_isothiocyanate_shape = _ISOTHIOCYANATE.has_shape
name_isothiocyanate = _ISOTHIOCYANATE.name
