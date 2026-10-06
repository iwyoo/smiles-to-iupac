"""'thiourea' (H2N-C(=S)-NH2) and its N-substituted derivatives, implemented in `_chalcogenourea.py`."""

from ._chalcogenourea import Chalcogenourea

_THIOUREA = Chalcogenourea(
    atomic_num=16,
    symbol="S",
    word="thiourea",
)

has_thiourea_shape = _THIOUREA.has_shape
name_thiourea = _THIOUREA.name
