"""Acyl prefixes of imidic and hydrazonic acids (P-65.1.7.2.2, P-65.1.7.4.2): ethanimidoyl, benzenecarboximidoyl,
N-hydroxybenzenecarboximidoyl, methanehydrazonoyl, formed from the name of the imidic acid made of the group."""

from rdkit import Chem

from ._common import UnsupportedStructure


def _subtree(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def imidoyl_prefix(mol, graph, root, coming_from):
    """(name, True) for a carbon group -C(R)=N-X, else None."""
    atom = mol.GetAtomWithIdx(root)
    if atom.IsInRing() or atom.GetFormalCharge() or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        return None
    others = [n for n in graph[root] if n != coming_from]
    imino = [
        n
        for n in others
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0
    ]
    rest = [n for n in others if n not in imino]
    if len(imino) != 1 or len(rest) > 1 or atom.GetTotalNumHs() != (0 if rest else 1):
        return None
    nitrogen = mol.GetAtomWithIdx(imino[0])
    if nitrogen.GetFormalCharge() or nitrogen.IsInRing():
        return None
    if rest and (
        mol.GetAtomWithIdx(rest[0]).GetAtomicNum() != 6
        or mol.GetBondBetweenAtoms(root, rest[0]).GetBondTypeAsDouble() != 1.0
    ):
        return None
    if any(mol.GetAtomWithIdx(n).GetAtomicNum() not in (1, 6, 7, 8) for n in graph[imino[0]] if n != root):
        return None
    if any(
        mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() > 1
        for n in graph[imino[0]]
        if n != root
    ):
        return None
    from ._functional_prefixes import _acyl_prefix

    subtree = _subtree(graph, root, coming_from)
    try:
        name = _acyl_prefix(mol, subtree, root, coming_from)
    except (UnsupportedStructure, Chem.rdchem.MolSanitizeException):
        return None
    return (name, True) if name.endswith(("imidoyl", "hydrazonoyl")) else None


_KETENE_PREFIX = {8: "oxo", 16: "sulfanylidene", 34: "selanylidene", 52: "tellanylidene"}


def ketene_prefix(mol, graph, root, coming_from):
    """(name, True) for the nonacyl group =C=X of carbonic acid analogues (P-65.2.1.8): oxomethylidene,
    sulfanylidenemethylidene."""
    atom = mol.GetAtomWithIdx(root)
    if atom.IsInRing() or atom.GetFormalCharge() or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 2.0:
        return None
    others = [n for n in graph[root] if n != coming_from]
    if len(others) != 1:
        return None
    end = mol.GetAtomWithIdx(others[0])
    if (
        end.GetAtomicNum() not in _KETENE_PREFIX
        or end.GetDegree() != 1
        or end.GetFormalCharge()
        or mol.GetBondBetweenAtoms(root, others[0]).GetBondTypeAsDouble() != 2.0
    ):
        return None
    return _KETENE_PREFIX[end.GetAtomicNum()] + "methylidene", True
