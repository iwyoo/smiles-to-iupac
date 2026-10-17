"""The ester of sulfuric, sulfurous, nitric or nitrous acid as the parent beside the groups that rank below an ester
(P-41, P-67.1.3.2): acid halides, amides, nitriles, aldehydes, ketones, alcohols and their chalcogen analogues, amines and
imines are prefixes of the organyl group."""

from ._common import UnsupportedStructure
from ._nitrate_ester import has_nitrate_ester_shape, name_nitrate_ester
from ._nitrite_ester import has_nitrite_ester_shape, name_nitrite_ester
from ._sulfate import has_sulfate_shape, name_sulfate
from ._sulfite import has_sulfite_shape, name_sulfite

_ESTERS = (
    (has_sulfate_shape, name_sulfate),
    (has_sulfite_shape, name_sulfite),
    (has_nitrate_ester_shape, name_nitrate_ester),
    (has_nitrite_ester_shape, name_nitrite_ester),
)


def name_inorganic_ester_parent(mol):
    for has_shape, namer in _ESTERS:
        if has_shape(mol):
            return namer(mol)
    raise UnsupportedStructure("no ester of sulfuric, sulfurous, nitric or nitrous acid")
