"""An amine of the lambda2 methane (P-61.10): R-NH-CH: is N-R-λ2-methanamine, N-hydroxy-λ2-methanamine for HO-NH-CH:."""

from ._common import adjacency, halogen_substituents
from ._substituents import name_branch


def _found(mol):
    if len(mol.GetAtoms()) < 3:
        return None
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetNumRadicalElectrons() == 2 and a.GetDegree() == 1 and a.GetTotalNumHs() == 1]
    if len(carbons) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    (carbon,) = carbons
    (nitrogen,) = carbon.GetNeighbors()
    if nitrogen.GetAtomicNum() != 7 or nitrogen.GetDegree() != 2 or nitrogen.GetTotalNumHs() != 1:
        return None
    (root,) = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
    if root.GetAtomicNum() not in (6, 8) or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), root.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    return nitrogen, root


def has_carbene_amine_shape(mol) -> bool:
    return _found(mol) is not None


def name_carbene_amine(mol) -> str:
    nitrogen, root = _found(mol)
    graph = adjacency(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    name, compound = name_branch(graph, root.GetIdx(), nitrogen.GetIdx(), halogen_substituents(mol), aromatic, mol=mol, unsaturated=True)
    return f"N-({name})-λ2-methanamine" if compound else f"N-{name}-λ2-methanamine"
