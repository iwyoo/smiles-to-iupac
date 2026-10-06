"""'R tellurocyanate' (R-Te-C#N), implemented in `_chalcocyanate.py`."""

from ._chalcocyanate import Cyanate

_TELLUROCYANATE = Cyanate(
    atomic_num=52,
    symbol="Te",
    element="tellurium",
    word="tellurocyanate",
)

has_tellurocyanate_shape = _TELLUROCYANATE.has_shape
name_tellurocyanate = _TELLUROCYANATE.name
