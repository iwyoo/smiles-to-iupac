"""'telluroxide' (P-63.6), implemented in `_chalcogen_oxide.py`."""

from ._chalcogen_oxide import ChalcogenOxide

_TELLUROXIDE = ChalcogenOxide(
    atomic_num=52,
    element="tellurium",
    word="telluroxide",
    prefix="tellurinyl",
    oxygens=1,
)

has_telluroxide_shape = _TELLUROXIDE.has_shape
name_telluroxide = _TELLUROXIDE.name
