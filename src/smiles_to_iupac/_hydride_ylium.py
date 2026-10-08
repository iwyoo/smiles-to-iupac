"""Cations of mononuclear hydrides formed by loss of a hydride ion (P-73.2.2.1): the parent hydride name with the ending
'ylium', its substituents cited in front: azanylium, phosphanylium, phenylsulfanylium, triphenylsilylium,
chloranylium. The cationic atom has one bond fewer than its standard bonding number."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_mononuclear_prefixes, name_branch

_STEM = {
    5: ("boran", 3), 7: ("azan", 3), 8: ("oxidan", 2), 14: ("sil", 4), 15: ("phosphan", 3), 16: ("sulfan", 2),
    17: ("chloran", 1), 32: ("germ", 4), 33: ("arsan", 3), 34: ("selan", 2), 35: ("broman", 1), 50: ("stann", 4),
    51: ("stiban", 3), 52: ("tellan", 2), 53: ("iodan", 1), 82: ("plumb", 4), 83: ("bismuthan", 3),
}


def _center(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or charged[0].GetFormalCharge() != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    atom = charged[0]
    if atom.GetAtomicNum() not in _STEM or atom.IsInRing() or atom.GetIsotope():
        return None
    standard = _STEM[atom.GetAtomicNum()][1]
    bonds = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()
    if bonds != standard - 1:
        return None
    if atom.GetAtomicNum() == 7 and atom.GetDegree():
        return None
    if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
        return None
    if any(n.GetAtomicNum() not in (6, 9, 17, 35, 53) for n in atom.GetNeighbors()):
        return None
    if any(a.GetIsotope() for a in mol.GetAtoms()) or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms() if a.GetIdx() != atom.GetIdx()):
        return None
    return atom


def has_hydride_ylium_shape(mol) -> bool:
    return _center(mol) is not None


def name_hydride_ylium(mol) -> str:
    center = _center(mol)
    if center is None:
        raise UnsupportedStructure("not the cation of a mononuclear hydride")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [
        name_branch(graph, n, center.GetIdx(), halogens, aromatic, mol=mol) for n in graph[center.GetIdx()]
    ]
    z = center.GetAtomicNum()
    stem = _STEM[z][0]
    parent = stem + "ylium"
    return format_mononuclear_prefixes(entries) + parent if entries else parent


_ONIUM_STEM = {7: "azanium", 8: "oxidanium", 16: "sulfanium", 34: "selanium", 52: "telluranium", 17: "chloranium", 35: "bromanium", 53: "iodanium"}


def _onium_center(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or charged[0].GetFormalCharge() != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    atom = charged[0]
    if atom.GetAtomicNum() not in _ONIUM_STEM or atom.IsInRing() or atom.GetIsotope():
        return None
    bonds = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()
    if bonds != _STEM[atom.GetAtomicNum()][1] + 1 or not atom.GetDegree():
        return None
    if any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors()) or any(b.GetBondTypeAsDouble() > 2.0 for b in atom.GetBonds()):
        return None
    if atom.GetAtomicNum() == 7 and all(b.GetBondTypeAsDouble() == 1.0 for b in atom.GetBonds()):
        return None
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms() if a.GetIdx() != atom.GetIdx()):
        return None
    return atom


def has_hydride_onium_shape(mol) -> bool:
    return _onium_center(mol) is not None


def name_hydride_onium(mol) -> str:
    """P-73.1.2.1: an oxonium, sulfonium or halonium centre with ylidene and alkyl groups: ethyl(propan-2-ylidene)oxidanium,
    acetyl(methyl)chloranium."""
    from ._dipolar import _group

    center = _onium_center(mol)
    if center is None:
        raise UnsupportedStructure("not an onium cation of a mononuclear hydride")
    graph = adjacency(mol)
    entries = []
    for n in graph[center.GetIdx()]:
        name, compound = _group(mol, graph, n, center.GetIdx())
        entries.append((name, compound))
    return format_mononuclear_prefixes(entries) + _ONIUM_STEM[center.GetAtomicNum()]
