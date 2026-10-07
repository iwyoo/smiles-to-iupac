"""The free mononuclear noncarbon oxoacids of P-67.1.1.1 (boric, nitric, phosphoric, sulfuric, perchloric acid, ...),
named from the central element, its oxo groups, its hydroxy groups and the hydrogens left on it; substituted acids are
named from these (P-67.1.1.2)."""

from rdkit import Chem

from ._common import UnsupportedStructure

# (atomic number, oxo, hydroxy, hydrogen) -> name
_ACIDS = {
    (5, 0, 3, 0): "boric acid",
    (5, 0, 2, 1): "boronic acid",
    (5, 0, 1, 2): "borinic acid",
    (7, 1, 3, 0): "nitroric acid",
    (7, 0, 3, 0): "azorous acid",
    (7, 1, 2, 1): "azonic acid",
    (7, 1, 1, 2): "azinic acid",
    (7, 0, 2, 1): "azonous acid",
    (7, 2, 1, 0): "nitric acid",
    (7, 1, 1, 0): "nitrous acid",
    (14, 0, 4, 0): "silicic acid",
    (16, 2, 2, 0): "sulfuric acid",
    (16, 1, 2, 0): "sulfurous acid",
    (34, 2, 2, 0): "selenic acid",
    (34, 1, 2, 0): "selenous acid",
    (52, 2, 2, 0): "telluric acid",
    (52, 1, 2, 0): "tellurous acid",
}
for _z, _stem in ((15, "phosphor"), (33, "arsor"), (51, "stibor")):
    _ACIDS[(_z, 1, 3, 0)] = f"{_stem}ic acid"
    _ACIDS[(_z, 0, 3, 0)] = f"{_stem}ous acid"
for _z, _stem in ((15, "phosphon"), (33, "arson"), (51, "stibon")):
    _ACIDS[(_z, 1, 2, 1)] = f"{_stem}ic acid"
    _ACIDS[(_z, 0, 2, 1)] = f"{_stem}ous acid"
for _z, _stem in ((15, "phosphin"), (33, "arsin"), (51, "stibin")):
    _ACIDS[(_z, 1, 1, 2)] = f"{_stem}ic acid"
    _ACIDS[(_z, 0, 1, 2)] = f"{_stem}ous acid"
for _z, _stem in ((9, "fluor"), (17, "chlor"), (35, "brom"), (53, "iod")):
    _ACIDS[(_z, 0, 1, 0)] = f"hypo{_stem}ous acid"
    _ACIDS[(_z, 1, 1, 0)] = f"{_stem}ous acid"
    _ACIDS[(_z, 2, 1, 0)] = f"{_stem}ic acid"
    _ACIDS[(_z, 3, 1, 0)] = f"per{_stem}ic acid"


def _counts(mol):
    """(central atomic number, oxo, hydroxy, hydrogen, anionic oxygens) or None when the molecule is not a free
    mononuclear acid or one of its oxoanions."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() != 8]
    if len(centers) != 1 or centers[0].GetIsotope():
        return None
    center = centers[0]
    oxo = hydroxy = ionic = 0
    for oxygen in center.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(center.GetIdx(), oxygen.GetIdx())
        if oxygen.GetDegree() != 1 or oxygen.GetIsotope():
            return None
        if bond.GetBondTypeAsDouble() == 2.0 and not oxygen.GetFormalCharge():
            oxo += 1
        elif bond.GetBondTypeAsDouble() == 1.0 and oxygen.GetFormalCharge() == -1 and center.GetFormalCharge() > 0:
            oxo += 1
        elif bond.GetBondTypeAsDouble() == 1.0 and oxygen.GetFormalCharge() == -1 and not oxygen.GetTotalNumHs():
            ionic += 1
        elif bond.GetBondTypeAsDouble() == 1.0 and not oxygen.GetFormalCharge() and oxygen.GetTotalNumHs() == 1:
            hydroxy += 1
        else:
            return None
    if mol.GetNumAtoms() != 1 + oxo + hydroxy + ionic:
        return None
    return center.GetAtomicNum(), oxo, hydroxy + ionic, center.GetTotalNumHs(), ionic


_IRREGULAR = {
    "sulfuric": "sulfate",
    "sulfurous": "sulfite",
    "phosphoric": "phosphate",
    "phosphorous": "phosphite",
}


def _anion_name(acid, hydrogens):
    stem = acid[: -len(" acid")]
    if stem in _IRREGULAR:
        stem = _IRREGULAR[stem]
    elif stem.endswith("ic"):
        stem = stem[:-2] + "ate"
    else:
        stem = stem[:-3] + "ite"
    words = {0: "", 1: "hydrogen ", 2: "dihydrogen "}[hydrogens]
    return words + stem


def has_mononuclear_oxoacid_shape(mol) -> bool:
    key = _counts(mol)
    return key is not None and key[:4] in _ACIDS and _balanced(mol, key)


def _balanced(mol, key):
    charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    return charge == -key[4]


def name_mononuclear_oxoacid(mol) -> str:
    key = _counts(mol)
    if key is None or key[:4] not in _ACIDS or not _balanced(mol, key):
        raise UnsupportedStructure("this is not a free mononuclear noncarbon oxoacid")
    acid = _ACIDS[key[:4]]
    if not key[4]:
        return acid
    return _anion_name(acid, key[2] - key[4])
