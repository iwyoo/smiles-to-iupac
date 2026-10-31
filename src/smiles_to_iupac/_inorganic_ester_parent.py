"""The ester of sulfuric, sulfurous, nitric, nitrous or a halogen oxoacid as the parent beside the groups that rank below an ester
(P-41, P-67.1.3.2): acid halides, amides, nitriles, aldehydes, ketones, alcohols and their chalcogen analogues, amines and
imines are prefixes of the organyl group."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._halogen_acid_ester import _ester_halogens, name_halogen_acid_ester
from ._nitrate_ester import has_nitrate_ester_shape, name_nitrate_ester
from ._nitrite_ester import has_nitrite_ester_shape, name_nitrite_ester
from ._sulfate import has_sulfate_shape, has_sulfite_shape, name_sulfate, name_sulfite

def _has_halogen_ester_parent(mol):
    """A cation or anion elsewhere in the molecule outranks the ester (P-41)."""
    groups = _ester_halogens(mol)
    if len(groups) != 1:
        return False
    halogen = groups[0][0]
    group = {halogen.GetIdx(), *(n.GetIdx() for n in halogen.GetNeighbors())}
    return not any(a.GetFormalCharge() and a.GetIdx() not in group for a in mol.GetAtoms())


def _cation_outranks_ester(mol):
    """A cationic centre in the same fragment, with no anion beside it, outranks the ester of a sulfur oxoacid (P-41)."""
    for fragment in Chem.GetMolFrags(mol):
        charges = [mol.GetAtomWithIdx(i).GetFormalCharge() for i in fragment]
        if any(c > 0 for c in charges) and not any(c < 0 for c in charges):
            return True
    return False


def _without_cation(has_shape):
    return lambda mol: has_shape(mol) and not _cation_outranks_ester(mol)


_ESTERS = (
    (_without_cation(has_sulfate_shape), name_sulfate),
    (_without_cation(has_sulfite_shape), name_sulfite),
    (has_nitrate_ester_shape, name_nitrate_ester),
    (has_nitrite_ester_shape, name_nitrite_ester),
    (_has_halogen_ester_parent, name_halogen_acid_ester),
)


def name_inorganic_ester_parent(mol):
    for has_shape, namer in _ESTERS:
        if has_shape(mol):
            return namer(mol)
    raise UnsupportedStructure("no ester of a sulfur, nitrogen or halogen oxoacid")


def has_inorganic_ester_parent(mol):
    return any(has_shape(mol) for has_shape, _ in _ESTERS)

