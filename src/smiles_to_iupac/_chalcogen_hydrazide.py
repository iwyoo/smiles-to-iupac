"""Chalcogen analogues of hydrazides (P-66.3.4): every hydrazide group is named with its =O in place of the =S, =Se or
=Te and the replacement infix 'thio', 'seleno' or 'telluro' is put in front of 'hydrazide'; the retained acetic,
benzoic and formic names are not used. Groups of one chalcogen are multiplied like the oxygen hydrazides."""

from rdkit import Chem

from ._common import UnsupportedStructure

_INFIX = {16: "thio", 34: "seleno", 52: "telluro"}
_HYDRAZIDE = Chem.MolFromSmarts("[CX3;$([CX3]-[#6]),$([CX3;H1]),$([CX3]-[NX3;R])](=[S,Se,Te;X1])-[NX3]-[NX3]")
_OXYGEN_HYDRAZIDE = Chem.MolFromSmarts("[CX3;$([CX3]-[#6]),$([CX3;H1]),$([CX3]-[NX3;R])](=[O;X1])-[NX3]-[NX3]")
_RETAINED = (("acetohydrazide", "ethane{}hydrazide"), ("benzohydrazide", "benzenecarbo{}hydrazide"), ("formohydrazide", "methane{}hydrazide"))


def has_chalcogen_hydrazide_shape(mol) -> bool:
    return mol.HasSubstructMatch(_HYDRAZIDE)


def name_chalcogen_hydrazide(mol) -> str:
    matches = mol.GetSubstructMatches(_HYDRAZIDE)
    chalcogens = {match[1] for match in matches}
    if len({mol.GetAtomWithIdx(a).GetAtomicNum() for a in chalcogens}) != 1:
        raise UnsupportedStructure("hydrazide groups of different chalcogens are not supported yet")
    if mol.HasSubstructMatch(_OXYGEN_HYDRAZIDE):
        raise UnsupportedStructure("a hydrazide group beside a chalcogen hydrazide group is named by the senior hydrazide")
    if any(a.GetIdx() not in chalcogens and a.GetAtomicNum() in _INFIX and a.GetDegree() == 1 for a in mol.GetAtoms()):
        raise UnsupportedStructure("further thione-type groups are not supported yet")
    infix = _INFIX[mol.GetAtomWithIdx(next(iter(chalcogens))).GetAtomicNum()]
    healed = Chem.RWMol(mol)
    for atom in chalcogens:
        healed.GetAtomWithIdx(atom).SetAtomicNum(8)
    healed = healed.GetMol()
    Chem.SanitizeMol(healed)
    from .core import _smiles_to_iupac_unabridged

    name = _smiles_to_iupac_unabridged(Chem.MolToSmiles(healed))
    if not name.endswith("hydrazide") or name.count("hydrazide") != 1:
        raise UnsupportedStructure("the oxygen analogue is not named by one hydrazide suffix")
    for retained, systematic in _RETAINED:
        if name.endswith(retained):
            return name[: -len(retained)] + systematic.format(infix)
    return name[: -len("hydrazide")] + infix + "hydrazide"
