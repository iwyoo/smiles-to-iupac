"""'selenoate' anions (P-72.2.2.2.1.1), implemented in `_chalcogenoate.py`."""

from ._chalcogenoate import Ate

_SELENOATE = Ate(
    atomic_num=34,
    symbol="Se",
    word="selenoate",
)

has_selenoate_shape = _SELENOATE.has_shape
name_selenoate = _SELENOATE.name
