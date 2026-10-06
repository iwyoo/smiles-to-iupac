"""'sulfinic acid' (P-65.3.1), implemented in `_oxo_acid.py`."""

from ._oxo_acid import OxoAcid

_SULFINIC_ACID = OxoAcid(
    atomic_num=16,
    word="sulfinic",
    element="sulfur",
    formula="SO2H",
    oxygens=2,
)

has_sulfinic_acid_shape = _SULFINIC_ACID.has_shape
name_sulfinic_acid = _SULFINIC_ACID.name
_sulfinic_sulfur_atoms = _SULFINIC_ACID._hetero_atoms
