"""'tellurone' (P-63.6), implemented in `_chalcogen_oxide.py`."""

from ._chalcogen_oxide import ChalcogenOxide

_TELLURONE = ChalcogenOxide(
    atomic_num=52,
    element="tellurium",
    word="tellurone",
    prefix="telluronyl",
    oxygens=2,
)

has_tellurone_shape = _TELLURONE.has_shape
name_tellurone = _TELLURONE.name
