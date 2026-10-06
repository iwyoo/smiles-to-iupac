"""'disulfide' and longer chalcogen chains R-S(n)-R', implemented in `_dichalcogenide.py`."""

from ._dichalcogenide import Dichalcogenide

_DISULFIDE = Dichalcogenide(
    atomic_num=16,
    symbol="S",
    element="sulfur",
    word="disulfide",
    base="sulfide",
    stem="sulfanyl",
    perol="perthiol",
)

has_disulfide_shape = _DISULFIDE.has_shape
name_disulfide = _DISULFIDE.name
