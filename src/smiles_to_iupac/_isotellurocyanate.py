"""'R-N=C=Te' (isotellurocyanate, P-61.8), implemented in `_iso_cyanate.py`."""

from ._iso_cyanate import Isocyanate

_ISOTELLUROCYANATE = Isocyanate(
    atomic_num=52,
    symbol="Te",
    element="tellurium",
    infix="telluro",
)

has_isotellurocyanate_shape = _ISOTELLUROCYANATE.has_shape
name_isotellurocyanate = _ISOTELLUROCYANATE.name
