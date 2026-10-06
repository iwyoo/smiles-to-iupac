"""'sulfonic acid' (P-65.3.1), implemented in `_oxo_acid.py`."""

from ._oxo_acid import OxoAcid

_SULFONIC_ACID = OxoAcid(
    atomic_num=16,
    word="sulfonic",
    element="sulfur",
    formula="SO3H",
    oxygens=3,
)

has_sulfonic_acid_shape = _SULFONIC_ACID.has_shape
name_sulfonic_acid = _SULFONIC_ACID.name
_sulfonic_sulfur_atoms = _SULFONIC_ACID._hetero_atoms
_name_acyclic_sulfonic_acid = _SULFONIC_ACID._name_acyclic
