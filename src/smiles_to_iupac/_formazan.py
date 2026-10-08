"""Formazan, HN=N-CH=N-NH2, a retained parent compound (P-52.2.2.1, P-68.3.1.3.5): the chain is numbered N1=N2-C3=N4-N5
from the azo end, so '1,3,5-triphenylformazan' and '1,5-diphenylformazan-3-carboxylic acid'. Substituents on the chain
are hydrocarbyl or halogen groups; a carboxylic acid on C3 is the suffix. Its substituent groups formazan-1-yl,
formazan-3-yl and formazan-5-yl, and the linker formazan-1,5-diyl and formazan-3,5-diyl, are the preferred prefixes
(P-34.2.1.3, P-68.3.1.3.5.2)."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents
from ._substituents import format_substituent_prefixes, name_branch

_CHAIN = Chem.MolFromSmarts("[NX2;!R]=[NX2;!R]-[CX3;!R]=[NX2;!R]-[NX3;!R]")
FORMAZAN_SKELETON = Chem.MolFromSmarts("[#7;!R]=[#7;!R]-[#6;!R]=[#7;!R]-[#7;!R]")
_GROUP_ATOMS = (6, 9, 17, 35, 53)


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


def _subtree(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n != blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _chain_groups(mol, graph, halogens, aromatic, chain, free):
    """{position: [(name, is_compound)]} of the hydrocarbyl and halogen groups on the chain atoms; None when another group
    is present. `free` holds the atoms outside the chain that are the free valences."""
    named = {}
    chain_set = set(chain)
    for position, atom in enumerate(chain, start=1):
        for n in graph[atom]:
            if n in chain_set or n in free:
                continue
            if mol.GetAtomWithIdx(n).GetAtomicNum() not in _GROUP_ATOMS or mol.GetBondBetweenAtoms(atom, n).GetBondTypeAsDouble() != 1.0:
                return None
            named.setdefault(position, []).append(name_branch(graph, n, atom, halogens, aromatic, mol=mol, unsaturated=True))
    return named


def formazan_substituent(mol, graph, root, coming_from, halogens, aromatic_atoms):
    """(name, True) of a formazan group joined through N1, C3 or N5, 'formazan-1-yl' or '(1,5-diphenylformazan-3-yl)'
    (P-34.2.1.3, P-68.3.1.3.5.2); None when the branch is anything else."""
    atoms = _subtree(graph, root, coming_from)
    if coming_from in atoms or any(mol.GetAtomWithIdx(a).GetFormalCharge() or mol.GetAtomWithIdx(a).GetIsotope() for a in atoms):
        return None
    for chain in mol.GetSubstructMatches(FORMAZAN_SKELETON):
        if root not in (chain[0], chain[2], chain[4]) or not set(chain) <= atoms:
            continue
        named = _chain_groups(mol, graph, halogens, aromatic_atoms, chain, {coming_from})
        if named is None:
            return None
        prefixes = format_substituent_prefixes(group_substituents(named)) if named else ""
        return f"{prefixes}formazan-{chain.index(root) + 1}-yl", True
    return None


def formazan_linker(mol, linker_atoms, units):
    """Text of the formazan skeleton joining `units` ([(junction atom of the unit, linker atom bonded to it)]) as
    'formazan-1,5-diyl', 'formazan-3,5-diyl' or 'formazan-1,3,5-triyl' with its substituent prefixes, or None."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    attached = {linker: junction for junction, linker in units}
    for chain in mol.GetSubstructMatches(FORMAZAN_SKELETON):
        if not set(chain) <= linker_atoms or any(a not in chain for a in attached):
            continue
        positions = sorted(chain.index(a) + 1 for a in attached)
        if positions not in ([1, 5], [3, 5], [1, 3, 5]) or len(attached) != len(units):
            continue
        if any(mol.GetBondBetweenAtoms(j, a).GetBondTypeAsDouble() != 1.0 for a, j in attached.items()):
            continue
        named = _chain_groups(mol, graph, halogens, aromatic, chain, set(attached.values()))
        if named is None:
            continue
        prefixes = format_substituent_prefixes(group_substituents(named)) if named else ""
        valence = {2: "diyl", 3: "triyl"}[len(positions)]
        return f"{prefixes}formazan-{','.join(map(str, positions))}-{valence}", bool(prefixes)
    return None
