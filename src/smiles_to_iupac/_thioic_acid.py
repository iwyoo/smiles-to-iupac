"""'thioic_acid' (P-65.1.5), implemented in `_chalcogenoic_acid.py`."""

from ._chalcogenoic_acid import ChalcogenoicAcid

_THIOIC_ACID = ChalcogenoicAcid(
    atomic_num=16,
    symbol="S",
    word="thioic",
)

has_thioic_acid_shape = _THIOIC_ACID.has_shape
name_thioic_acid = _THIOIC_ACID.name
