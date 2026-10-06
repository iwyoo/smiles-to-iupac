"""A carboxylic acid combined with sulfinic acid groups, implemented in `_carboxylic_acid_oxo_acid.py`."""

from ._carboxylic_acid_oxo_acid import CarboxylicAcidCombo
from ._sulfinic_acid import _SULFINIC_ACID

_CARBOXYLIC_ACID_SULFINIC_ACID = CarboxylicAcidCombo(
    atomic_num=16,
    word="sulfinic",
    prefix="sulfino",
    element="sulfur",
    acid=_SULFINIC_ACID,
)

has_carboxylic_acid_sulfinic_acid_shape = _CARBOXYLIC_ACID_SULFINIC_ACID.has_shape
name_carboxylic_acid_sulfinic_acid = _CARBOXYLIC_ACID_SULFINIC_ACID.name
