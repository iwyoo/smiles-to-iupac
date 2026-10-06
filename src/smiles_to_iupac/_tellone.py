"""'-tellone' (C=Te), the P-64.6.1 chalcogen analogue of a ketone; implemented in `_chalcogenone.py`."""

from ._chalcogenone import Chalcogenone

_TELLONE = Chalcogenone(
    atomic_num=52,
    symbol="Te",
    element="tellurium",
    suffix="tellone",
    ylidene="tellanylidene",
    ether="telluroether/telluride",
    hydride="tellurol, -TeH",
    aldehyde="telluroaldehyde",
)

has_tellone_shape = _TELLONE.has_shape
name_tellone = _TELLONE.name
