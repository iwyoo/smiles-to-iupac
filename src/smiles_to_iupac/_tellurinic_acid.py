"""'tellurinic acid' (P-65.3.1), implemented in `_oxo_acid.py`."""

from ._oxo_acid import OxoAcid

_TELLURINIC_ACID = OxoAcid(
    atomic_num=52,
    word="tellurinic",
    element="tellurium",
    formula="Te(=O)OH",
    oxygens=2,
)

has_tellurinic_acid_shape = _TELLURINIC_ACID.has_shape
name_tellurinic_acid = _TELLURINIC_ACID.name
