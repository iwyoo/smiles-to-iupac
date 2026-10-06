"""'telluronic acid' (P-65.3.1), implemented in `_oxo_acid.py`."""

from ._oxo_acid import OxoAcid

_TELLURONIC_ACID = OxoAcid(
    atomic_num=52,
    word="telluronic",
    element="tellurium",
    formula="Te(=O)(=O)OH",
    oxygens=3,
)

has_telluronic_acid_shape = _TELLURONIC_ACID.has_shape
name_telluronic_acid = _TELLURONIC_ACID.name
