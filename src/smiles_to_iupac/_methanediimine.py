"""Carbodiimides named as methanediimines (P-62.3.1.4): R-N=C=N-R' is the N,N'-substituted derivative of methanediimine,
the substituents of two identical nitrogen groups being cited without locants like those of any mononuclear hydride."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, group_substituents, halogen_substituents
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch


def _nitrogens(mol):
    """The two imino nitrogens of the single C(=N)=N carbon, else None."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    centers = [
        a for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetDegree() == 2 and not a.IsInRing()
        and all(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(a).GetAtomicNum() == 7 for b in a.GetBonds())
    ]
    if len(centers) != 1:
        return None
    nitrogens = [n for n in centers[0].GetNeighbors()]
    if any(n.GetDegree() > 2 or n.GetFormalCharge() or n.IsInRing() for n in nitrogens):
        return None
    return centers[0], nitrogens


def has_methanediimine_shape(mol) -> bool:
    found = _nitrogens(mol)
    if found is None:
        return False
    nitrogens = {n.GetIdx() for n in found[1]}
    return all(
        a.GetIdx() in nitrogens or (a.GetAtomicNum() in (6, *HALOGEN_PREFIXES) and not a.GetFormalCharge() and not a.GetIsotope())
        for a in mol.GetAtoms()
    )


def name_methanediimine(mol) -> str:
    center, nitrogens = _nitrogens(mol)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    cited = []
    for nitrogen in nitrogens:
        roots = [r for r in graph[nitrogen.GetIdx()] if r != center.GetIdx()]
        if any(mol.GetBondBetweenAtoms(nitrogen.GetIdx(), r).GetBondTypeAsDouble() != 1.0 for r in roots):
            raise UnsupportedStructure("a multiple bond on an imino nitrogen is not supported yet")
        cited.append([name_branch(graph, r, nitrogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True) for r in roots])
    if cited[0] == cited[1]:
        return format_mononuclear_prefixes(cited[0] + cited[1]) + "methanediimine"
    grouped = group_substituents(
        {locant: names for locant, names in (("N", cited[0]), ("N'", cited[1])) if names}
    )
    return format_substituent_prefixes(grouped) + "methanediimine"
