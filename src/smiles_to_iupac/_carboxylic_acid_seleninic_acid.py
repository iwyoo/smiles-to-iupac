"""A carboxylic acid combined with seleninic acid groups, implemented in `_carboxylic_acid_oxo_acid.py`."""

from ._carboxylic_acid_oxo_acid import CarboxylicAcidCombo
from ._seleninic_acid import _SELENINIC_ACID

_CARBOXYLIC_ACID_SELENINIC_ACID = CarboxylicAcidCombo(
    atomic_num=34,
    word="seleninic",
    prefix="selenino",
    element="selenium",
    acid=_SELENINIC_ACID,
)

has_carboxylic_acid_seleninic_acid_shape = _CARBOXYLIC_ACID_SELENINIC_ACID.has_shape
name_carboxylic_acid_seleninic_acid = _CARBOXYLIC_ACID_SELENINIC_ACID.name
