"""'telluroic_acid' (P-65.1.5), implemented in `_chalcogenoic_acid.py`."""

from ._chalcogenoic_acid import ChalcogenoicAcid

_TELLUROIC_ACID = ChalcogenoicAcid(
    atomic_num=52,
    symbol="Te",
    word="telluroic",
)

has_telluroic_acid_shape = _TELLUROIC_ACID.has_shape
name_telluroic_acid = _TELLUROIC_ACID.name
