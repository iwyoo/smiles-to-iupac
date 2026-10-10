"""Chalcogen analogues of hydrazides (P-66.3.4): the hydrazide is named with its =O in place of the =S, =Se or =Te and
the replacement infix 'thio', 'seleno' or 'telluro' is put in front of 'hydrazide'; the retained acetic, benzoic and
formic names are not used."""

from rdkit import Chem

from ._common import UnsupportedStructure

_INFIX = {16: "thio", 34: "seleno", 52: "telluro"}
_HYDRAZIDE = Chem.MolFromSmarts("[CX3;$([CX3]-[#6]),$([CX3;H1])](=[S,Se,Te;X1])-[NX3]-[NX3]")
_RETAINED = (("acetohydrazide", "ethane{}hydrazide"), ("benzohydrazide", "benzenecarbo{}hydrazide"), ("formohydrazide", "methane{}hydrazide"))


def has_chalcogen_hydrazide_shape(mol) -> bool:
    return mol.HasSubstructMatch(_HYDRAZIDE)


def name_chalcogen_hydrazide(mol) -> str:
    matches = mol.GetSubstructMatches(_HYDRAZIDE)
    if len(matches) != 1:
        raise UnsupportedStructure("exactly one chalcogen hydrazide group is required")
    carbon, chalcogen, _, _ = matches[0]
    if any(
        a.GetIdx() != chalcogen and a.GetAtomicNum() in _INFIX and a.GetDegree() == 1 for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("further thione-type groups are not supported yet")
    infix = _INFIX[mol.GetAtomWithIdx(chalcogen).GetAtomicNum()]
    healed = Chem.RWMol(mol)
    healed.GetAtomWithIdx(chalcogen).SetAtomicNum(8)
    healed = healed.GetMol()
    Chem.SanitizeMol(healed)
    from .core import _smiles_to_iupac_unabridged

    name = _smiles_to_iupac_unabridged(Chem.MolToSmiles(healed))
    if not name.endswith("hydrazide") or name.endswith("dihydrazide"):
        raise UnsupportedStructure("the oxygen analogue is not a single hydrazide")
    for retained, systematic in _RETAINED:
        if name.endswith(retained):
            return name[: -len(retained)] + systematic.format(infix)
    return name[: -len("hydrazide")] + infix + "hydrazide"
