"""'selenol' (R-SeH, P-63.1.1), implemented in `_chalcogenol.py`."""

from ._chalcogenol import Chalcogenol

_SELENOL = Chalcogenol(
    atomic_num=34,
    symbol="Se",
    element="selenium",
    word="selenol",
    phenol="selenophenol",
    base="selenide",
    stem="selanyl",
    article="an",
    coexisting="-OH, -SH, or amine",
)

has_selenol_shape = _SELENOL.has_shape
name_selenol = _SELENOL.name
