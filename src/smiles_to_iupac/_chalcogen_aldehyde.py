"""Chalcogen analogues of aldehydes (P-66.6.3, Table 6.4): the thial, selenal and tellanal groups are named by
naming the aldehyde and replacing its oxygen; the one-carbon retained names of P-66.5.1.2.1 and P-66.6.1.2.1 are
added."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure

# element -> (chain suffix, ring suffix tail after 'carbo', formyl-type prefix)
_GROUP = {
    16: ("thial", "thialdehyde", "methanethioyl"),
    34: ("selenal", "selenaldehyde", "methaneselenoyl"),
    52: ("tellanal", "telluraldehyde", "methanetelluroyl"),
}
_MULTIPLIERS = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
_WHOLE_MOLECULES = {"C=O": "formaldehyde", "C#N": "formonitrile", "C=S": "methanethial", "C=[Se]": "methaneselenal", "C=[Te]": "methanetellanal"}
_ALDEHYDE = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
_FORMYL_ACID = re.compile(r"(\d+)-formyl([a-z]+(?:-\d+-[a-z]+)? acid)")


def _groups(mol):
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() or atom.GetIsAromatic() or atom.GetTotalNumHs() != 1 or atom.GetDegree() != 2:
            continue
        for n in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx())
            if n.GetAtomicNum() in _GROUP and n.GetDegree() == 1 and not n.GetFormalCharge() and bond.GetBondTypeAsDouble() == 2.0:
                found.append(n.GetIdx())
    return found


def has_chalcogen_aldehyde_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _WHOLE_MOLECULES or bool(_groups(mol))


def _rename(name, element, count):
    chain, ring, prefix = _GROUP[element]
    multiplier = _MULTIPLIERS.get(count)
    if multiplier is None:
        raise UnsupportedStructure("too many chalcogen aldehyde groups")
    retained = name.replace("acetaldehyde", "ethane" + chain).replace("benzaldehyde", "benzenecarbo" + ring)
    if retained != name:
        return retained
    if name.endswith("carbaldehyde"):
        return name[: -len("carbaldehyde")] + "carbo" + ring
    ending = multiplier + "al"
    if name.endswith(ending):
        base = name[: -len(ending)]
        if not multiplier and base[-2:] in ("an", "en", "yn"):
            base += "e"
        return base + multiplier + chain
    match = _FORMYL_ACID.fullmatch(name)
    if match and count == 1:
        return f"{match.group(1)}-({prefix}){match.group(2)}"
    raise UnsupportedStructure("this chalcogen aldehyde is expressed by a prefix that is not supported yet")


def name_chalcogen_aldehyde(mol) -> str:
    whole = _WHOLE_MOLECULES.get(Chem.MolToSmiles(mol))
    if whole is not None:
        return whole
    chalcogens = _groups(mol)
    if not chalcogens:
        raise UnsupportedStructure("no chalcogen aldehyde group")
    if mol.HasSubstructMatch(_ALDEHYDE):
        raise UnsupportedStructure("mixed aldehyde and chalcogen aldehyde groups are not supported yet")
    senior = min((mol.GetAtomWithIdx(x).GetAtomicNum() for x in chalcogens))
    chalcogens = [x for x in chalcogens if mol.GetAtomWithIdx(x).GetAtomicNum() == senior]
    elements = {senior}
    healed = Chem.RWMol(mol)
    for x in chalcogens:
        healed.GetAtomWithIdx(x).SetAtomicNum(8)
    healed = healed.GetMol()
    Chem.SanitizeMol(healed)
    from .core import _smiles_to_iupac_unabridged

    healed_name = _smiles_to_iupac_unabridged(Chem.MolToSmiles(healed))
    original_oxo = any(
        a.GetAtomicNum() == 8 and a.GetDegree() == 1 and mol.GetBondBetweenAtoms(a.GetIdx(), a.GetNeighbors()[0].GetIdx()).GetBondTypeAsDouble() == 2.0
        for a in mol.GetAtoms()
    )
    if "oxo" in healed_name and not original_oxo:
        raise UnsupportedStructure("a chalcogen aldehyde group expressed as a prefix beside another one is not supported yet")
    return _rename(healed_name, elements.pop(), len(chalcogens))
