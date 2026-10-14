"""'selenourea' (H2N-C(=Se)-NH2) and its N-substituted derivatives, implemented in `_chalcogenourea.py`."""

from ._chalcogenourea import Chalcogenourea

_SELENOUREA = Chalcogenourea(atomic_num=34)

has_selenourea_shape = _SELENOUREA.has_shape
name_selenourea = _SELENOUREA.name
