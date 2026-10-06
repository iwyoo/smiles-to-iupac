"""'thioate' anions (P-72.2.2.2.1.1), implemented in `_chalcogenoate.py`."""

from ._chalcogenoate import Ate

_THIOATE = Ate(
    atomic_num=16,
    symbol="S",
    word="thioate",
)

has_thioate_shape = _THIOATE.has_shape
name_thioate = _THIOATE.name
