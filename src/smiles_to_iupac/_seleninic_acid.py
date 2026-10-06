"""'seleninic acid' (P-65.3.1), implemented in `_oxo_acid.py`."""

from ._oxo_acid import OxoAcid

_SELENINIC_ACID = OxoAcid(
    atomic_num=34,
    word="seleninic",
    element="selenium",
    formula="Se(=O)OH",
    oxygens=2,
)

has_seleninic_acid_shape = _SELENINIC_ACID.has_shape
name_seleninic_acid = _SELENINIC_ACID.name
_seleninic_selenium_atoms = _SELENINIC_ACID._hetero_atoms
