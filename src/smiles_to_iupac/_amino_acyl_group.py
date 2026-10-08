"""Acyl groups of the common amino acids as substituent prefixes: L-alanyl, glycyl, ... (P-103.2.5).

The group is closed with a hydroxy group, named as the free amino acid, and the ending is changed to 'yl'.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._peptide import acyl_word


def _is_amino_acyl_root(mol, graph, root, coming_from):
    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetFormalCharge() or atom.GetDegree() != 3:
        return False
    others = [n for n in graph[root] if n != coming_from]
    oxo = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0]
    if len(oxo) != 1:
        return False
    (alpha,) = [n for n in others if n != oxo[0]]
    alpha_atom = mol.GetAtomWithIdx(alpha)
    return alpha_atom.GetAtomicNum() == 6 and any(mol.GetAtomWithIdx(n).GetAtomicNum() == 7 for n in graph[alpha])


def amino_acyl_group(mol, graph, root, coming_from):
    """(name, True) for the acyl group at `root` when it is that of a common amino acid, else None."""
    from .core import smiles_to_iupac

    if not _is_amino_acyl_root(mol, graph, root, coming_from):
        return None
    inside, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != coming_from and n not in inside:
                inside.add(n)
                stack.append(n)
    if any(coming_from in graph[a] for a in inside if a != root):
        return None
    editable = Chem.RWMol(mol)
    oxygen = editable.AddAtom(Chem.Atom(8))
    editable.AddBond(root, oxygen, Chem.BondType.SINGLE)
    for index in sorted(set(range(mol.GetNumAtoms())) - inside, reverse=True):
        editable.RemoveAtom(index)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    acid = editable.GetMol()
    try:
        Chem.SanitizeMol(acid)
        name = acyl_word(smiles_to_iupac(Chem.MolToSmiles(acid)), xi=False)
    except (UnsupportedStructure, ValueError):
        return None
    return (name, True) if name else None
