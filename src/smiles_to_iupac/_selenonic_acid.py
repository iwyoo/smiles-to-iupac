"""'selenonic acid' (P-65.3.1), implemented in `_oxo_acid.py`."""

from ._oxo_acid import OxoAcid

_SELENONIC_ACID = OxoAcid(
    atomic_num=34,
    word="selenonic",
    element="selenium",
    formula="Se(=O)(=O)OH",
    oxygens=3,
)

has_selenonic_acid_shape = _SELENONIC_ACID.has_shape
name_selenonic_acid = _SELENONIC_ACID.name
