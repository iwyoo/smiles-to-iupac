"""A carboxylic acid combined with sulfonic acid groups, implemented in `_carboxylic_acid_oxo_acid.py`."""

from ._carboxylic_acid_oxo_acid import CarboxylicAcidCombo
from ._sulfonic_acid import _SULFONIC_ACID

_CARBOXYLIC_ACID_SULFONIC_ACID = CarboxylicAcidCombo(
    atomic_num=16,
    word="sulfonic",
    prefix="sulfo",
    element="sulfur",
    acid=_SULFONIC_ACID,
)

has_carboxylic_acid_sulfonic_acid_shape = _CARBOXYLIC_ACID_SULFONIC_ACID.has_shape
name_carboxylic_acid_sulfonic_acid = _CARBOXYLIC_ACID_SULFONIC_ACID.name
