"""P-21.1.1.2: water, ammonia and the binary names of the Group 16 and 17
hydrides are used, but the Blue Book defers their PIN status."""

from ._pin import mark

_COMMON_NAMES = {
    (7, 3): "ammonia",
    (8, 2): "water",
    (9, 1): "hydrogen fluoride",
    (16, 2): "hydrogen sulfide",
    (17, 1): "hydrogen chloride",
    (34, 2): "hydrogen selenide",
    (35, 1): "hydrogen bromide",
    (52, 2): "hydrogen telluride",
    (53, 1): "hydrogen iodide",
}

_DEFERRED = "the Blue Book defers the preferred name of this common-named hydride (P-21.1.1.2, P-14.8)"


def _common_name(mol):
    if mol.GetNumAtoms() != 1:
        return None
    (atom,) = mol.GetAtoms()
    if atom.GetFormalCharge() or atom.GetIsotope() or atom.GetNumRadicalElectrons():
        return None
    return _COMMON_NAMES.get((atom.GetAtomicNum(), atom.GetTotalNumHs()))


def has_common_hydride_shape(mol) -> bool:
    return _common_name(mol) is not None


def name_common_hydride(mol) -> str:
    return mark(_common_name(mol), _DEFERRED)
