"""'-thione' (C=S), the P-64.6.1 chalcogen analogue of a ketone; implemented in `_chalcogenone.py`."""

from ._chalcogenone import Chalcogenone

_THIONE = Chalcogenone(
    atomic_num=16,
    symbol="S",
    element="sulfur",
    suffix="thione",
    ylidene="sulfanylidene",
    ether="thioether/sulfide",
    hydride="thiol, -SH",
    aldehyde="thial/thioaldehyde",
)

has_thione_shape = _THIONE.has_shape
name_thione = _THIONE.name
