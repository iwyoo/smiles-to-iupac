"""'tellurol' (R-TeH, P-63.1.1), implemented in `_chalcogenol.py`."""

from ._chalcogenol import Chalcogenol

_TELLUROL = Chalcogenol(
    atomic_num=52,
    symbol="Te",
    element="tellurium",
    word="tellurol",
    phenol="tellurophenol",
    base="telluride",
    stem="tellanyl",
    article="a",
    coexisting="-OH, -SH, -SeH, or amine",
)

has_tellurol_shape = _TELLUROL.has_shape
name_tellurol = _TELLUROL.name
