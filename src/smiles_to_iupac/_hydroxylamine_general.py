"""O-substituted hydroxylamines with any carbon group on the oxygen (P-68.3.1.1.1.2): H2N-O-R is
'O-(chloromethyl)hydroxylamine', 'O-(cyclohexylmethyl)hydroxylamine', the group R named by the substituent machinery."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_substituent_prefixes, name_branch


def _atoms(mol):
    nitrogens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    if len(nitrogens) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    nitrogen = nitrogens[0]
    if nitrogen.GetDegree() != 1 or nitrogen.GetTotalNumHs() != 2 or nitrogen.GetFormalCharge():
        return None
    (oxygen,) = nitrogen.GetNeighbors()
    if oxygen.GetAtomicNum() != 8 or oxygen.GetDegree() != 2 or oxygen.GetFormalCharge() or oxygen.IsInRing():
        return None
    root = next(n for n in oxygen.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx())
    if root.GetAtomicNum() != 6:
        return None
    if any(a.GetAtomicNum() not in (6, 7, 8, 9, 17, 35, 53) or a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    if sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8) != 1:
        return None
    return nitrogen, oxygen, root


def has_o_substituted_hydroxylamine_shape(mol) -> bool:
    return _atoms(mol) is not None


def name_o_substituted_hydroxylamine(mol) -> str:
    found = _atoms(mol)
    if found is None:
        raise UnsupportedStructure("not an O-substituted hydroxylamine")
    nitrogen, oxygen, root = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    name, compound = name_branch(graph, root.GetIdx(), oxygen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
    return format_substituent_prefixes({name: {"locants": ["O"], "compound": compound}}) + "hydroxylamine"
