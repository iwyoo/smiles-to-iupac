"""'selenoic_acid' (P-65.1.5), implemented in `_chalcogenoic_acid.py`."""

from ._chalcogenoic_acid import ChalcogenoicAcid

_SELENOIC_ACID = ChalcogenoicAcid(
    atomic_num=34,
    symbol="Se",
    word="selenoic",
)

has_selenoic_acid_shape = _SELENOIC_ACID.has_shape
name_selenoic_acid = _SELENOIC_ACID.name
