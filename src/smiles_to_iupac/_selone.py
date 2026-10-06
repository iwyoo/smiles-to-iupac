"""'-selone' (C=Se), the P-64.6.1 chalcogen analogue of a ketone; implemented in `_chalcogenone.py`."""

from ._chalcogenone import Chalcogenone

_SELONE = Chalcogenone(
    atomic_num=34,
    symbol="Se",
    element="selenium",
    suffix="selone",
    ylidene="selanylidene",
    ether="selenoether/selenide",
    hydride="selenol, -SeH",
    aldehyde="selenoaldehyde",
)

has_selone_shape = _SELONE.has_shape
name_selone = _SELONE.name
