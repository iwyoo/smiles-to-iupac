"""Hydroxylamine and its chalcogen analogues H2N-E-R (E = O, S, Se, Te): O-substituted hydroxylamines with any carbon group
on the oxygen (P-68.3.1.1.1.2), 'O-(chloromethyl)hydroxylamine'; the O-carboxylic, O-carboxamide and O-sulfonic forms
(P-68.3.1.1.1.4), 'hydroxylamine-O-carboxylic acid'; and the thio, seleno and telluro analogues (P-68.3.1.1.1.6),
'S-methyl(thiohydroxylamine)'. The group R is named by the substituent machinery."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_substituent_prefixes, name_branch

_CHALCOGEN = {8: ("O", ""), 16: ("S", "thio"), 34: ("Se", "seleno"), 52: ("Te", "telluro")}
_ACID_FORMS = {
    "NOC(=O)O": "carboxylic acid",
    "NOC(N)=O": "carboxamide",
    "NOS(=O)(=O)O": "sulfonic acid",
}


def _acid_form(mol):
    return _ACID_FORMS.get(Chem.MolToSmiles(mol)) or _ACID_FORMS.get(Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol))))


def _atoms(mol):
    nitrogens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    if len(nitrogens) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    nitrogen = nitrogens[0]
    if nitrogen.GetDegree() != 1 or nitrogen.GetTotalNumHs() != 2 or nitrogen.GetFormalCharge():
        return None
    (chalcogen,) = nitrogen.GetNeighbors()
    if chalcogen.GetAtomicNum() not in _CHALCOGEN or chalcogen.GetFormalCharge() or chalcogen.IsInRing() or chalcogen.GetIsotope():
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    if chalcogen.GetDegree() == 1:
        return (nitrogen, chalcogen, None) if chalcogen.GetAtomicNum() != 8 else None
    if chalcogen.GetDegree() != 2:
        return None
    root = next(n for n in chalcogen.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx())
    if root.GetAtomicNum() != 6:
        return None
    allowed = (6, 7, 9, 17, 35, 53, chalcogen.GetAtomicNum())
    if any(a.GetAtomicNum() not in allowed for a in mol.GetAtoms()):
        return None
    if sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in _CHALCOGEN) != 1:
        return None
    return nitrogen, chalcogen, root


def has_o_substituted_hydroxylamine_shape(mol) -> bool:
    return _atoms(mol) is not None or _acid_form(mol) is not None


def name_o_substituted_hydroxylamine(mol) -> str:
    acid = _acid_form(mol)
    if acid is not None:
        return f"hydroxylamine-O-{acid}"
    found = _atoms(mol)
    if found is None:
        raise UnsupportedStructure("not an O-substituted hydroxylamine or an analogue")
    nitrogen, chalcogen, root = found
    symbol, word = _CHALCOGEN[chalcogen.GetAtomicNum()]
    parent = f"{word}hydroxylamine"
    if root is None:
        return parent
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    name, compound = name_branch(graph, root.GetIdx(), chalcogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
    prefixes = format_substituent_prefixes({name: {"locants": [symbol], "compound": compound}})
    return f"{prefixes}({parent})" if word else f"{prefixes}{parent}"
