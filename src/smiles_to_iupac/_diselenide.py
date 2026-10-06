"""'diselenide' (R-Se-Se-R'), implemented in `_dichalcogenide.py`."""

from ._dichalcogenide import Dichalcogenide

_DISELENIDE = Dichalcogenide(
    atomic_num=34,
    symbol="Se",
    element="selenium",
    word="diselenide",
    prefix="diselanyl",
)

has_diselenide_shape = _DISELENIDE.has_shape
name_diselenide = _DISELENIDE.name
