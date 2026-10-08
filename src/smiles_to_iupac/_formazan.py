"""Formazan, HN=N-CH=N-NH2, a retained parent compound (P-52.2.2.1, P-68.3.1.3.5): the chain is numbered N1=N2-C3=N4-N5
from the azo end, so '1,3,5-triphenylformazan' and '1,5-diphenylformazan-3-carboxylic acid'. Substituents on the chain
are hydrocarbyl or halogen groups; a carboxylic acid on C3 is the suffix."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents
from ._substituents import format_substituent_prefixes, name_branch

_CHAIN = Chem.MolFromSmarts("[NX2;!R]=[NX2;!R]-[CX3;!R]=[NX2;!R]-[NX3;!R]")


def _match(mol):
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    hits = mol.GetSubstructMatches(_CHAIN)
    if len(hits) != 1:
        return None
    n1, n2, c3, n4, n5 = hits[0]
    chain = set(hits[0])
    if any(a.GetAtomicNum() == 7 and a.GetIdx() not in chain for a in mol.GetAtoms()):
        return None
    return hits[0]


def has_formazan_shape(mol) -> bool:
    return _match(mol) is not None


def name_formazan(mol) -> str:
    chain = _match(mol)
    if chain is None:
        raise UnsupportedStructure("not a formazan")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    chain_set = set(chain)
    acid = None
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        for n in graph[atom]:
            if n in chain_set:
                continue
            neighbour = mol.GetAtomWithIdx(n)
            if position == 3 and neighbour.GetAtomicNum() == 6 and _is_carboxylic_acid(mol, n):
                if acid is not None:
                    raise UnsupportedStructure("several acid groups on a formazan are not supported yet")
                acid = n
                continue
            if neighbour.GetAtomicNum() not in (6, 9, 17, 35, 53) or mol.GetBondBetweenAtoms(atom, n).GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure("this substituent of formazan is not supported yet")
            substituents.setdefault(position, []).append(name_branch(graph, n, atom, halogens, aromatic, mol=mol, unsaturated=True))
    if acid is None and any(
        a.GetIdx() not in chain_set and a.GetAtomicNum() == 8 for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("an oxygen group on a formazan is not supported yet")
    grouped = group_substituents(substituents)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    name = f"{prefix}formazan"
    return f"{name}-3-carboxylic acid" if acid is not None else name


def _is_carboxylic_acid(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
    return (
        atom.GetDegree() == 3
        and len(oxygens) == 2
        and sorted(mol.GetBondBetweenAtoms(idx, o.GetIdx()).GetBondTypeAsDouble() for o in oxygens) == [1.0, 2.0]
        and any(o.GetTotalNumHs() == 1 for o in oxygens)
    )
