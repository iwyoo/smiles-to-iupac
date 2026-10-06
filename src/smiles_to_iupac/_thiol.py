"""'thiol' (R-SH, P-63.1.1), implemented in `_chalcogenol.py`."""

from ._chalcogenol import Chalcogenol

_THIOL = Chalcogenol(
    atomic_num=16,
    symbol="S",
    element="sulfur",
    word="thiol",
    phenol="thiophenol",
    base="sulfide",
    stem="sulfanyl",
    article="an",
    coexisting="-OH or amine",
)

has_thiol_shape = _THIOL.has_shape
name_thiol = _THIOL.name

_name_acyclic_thiol = _THIOL._name_acyclic
