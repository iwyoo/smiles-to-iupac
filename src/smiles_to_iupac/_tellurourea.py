"""'tellurourea' (H2N-C(=Te)-NH2) and its N-substituted derivatives, implemented in `_chalcogenourea.py`."""

from ._chalcogenourea import Chalcogenourea

_TELLUROUREA = Chalcogenourea(
    atomic_num=52,
    symbol="Te",
    word="tellurourea",
)

has_tellurourea_shape = _TELLUROUREA.has_shape
name_tellurourea = _TELLUROUREA.name
