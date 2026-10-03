"""Ester names built from their two parts (P-65.6.3.2.1): the alkyl/aryl group
cited as a substituent group, then the anion of the acid part. The molecule is
split at the O-alkyl bond; the acid fragment is named as a free acid and its
'-ic acid' ending becomes '-ate', so every acid shape the engine can name also
works as an ester acyl part.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import name_branch

_ESTER = Chem.MolFromSmarts("[CX3;!R](=O)[OX2;!R][#6]")
_FREE_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
_HYDROGEN_WORDS = {1: "hydrogen", 2: "dihydrogen", 3: "trihydrogen"}


def _carbonyl_with_heteroatom(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6:
        return False
    has_double_o = any(
        n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in atom.GetNeighbors()
    )
    return has_double_o and any(
        n.GetAtomicNum() in (7, 8, 16, 9, 17, 35, 53) and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in atom.GetNeighbors()
    )


def _anion_name(acid_name):
    if not acid_name.endswith(" acid") or " " in acid_name[: -len(" acid")].strip():
        raise UnsupportedStructure("the acid part of this ester is not named as a plain acid")
    stem = acid_name[: -len(" acid")]
    if not stem.endswith("ic"):
        raise UnsupportedStructure("the acid part of this ester has an unexpected acid name ending")
    return stem[:-2] + "ate"


def _branch_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def name_ester_by_parts(mol) -> str:
    from .core import smiles_to_iupac

    matches = mol.GetSubstructMatches(_ESTER)
    if len(matches) != 1:
        raise UnsupportedStructure("exactly one acyclic ester group is required for part-wise ester naming")
    acyl_carbon, _, ester_oxygen, alkyl_carbon = matches[0]
    graph = adjacency(mol)
    alkyl_atoms = _branch_atoms(graph, alkyl_carbon, ester_oxygen)
    if acyl_carbon in alkyl_atoms:
        raise UnsupportedStructure("a ring-closing ester (lactone) is not named part-wise")
    if any(_carbonyl_with_heteroatom(mol, a) for a in alkyl_atoms):
        raise UnsupportedStructure("an acid or ester group inside the alkyl part outranks this ester")
    if any(mol.GetAtomWithIdx(a).GetIsotope() or mol.GetAtomWithIdx(a).GetNumRadicalElectrons() for a in alkyl_atoms):
        raise UnsupportedStructure("isotopes and radicals in the alkyl part are not supported")

    editable = Chem.RWMol(mol)
    editable.RemoveBond(ester_oxygen, alkyl_carbon)
    oxygen = editable.GetAtomWithIdx(ester_oxygen)
    oxygen.SetNumExplicitHs(1)
    oxygen.SetNoImplicit(True)
    fragments = Chem.GetMolFrags(editable, asMols=True, sanitizeFrags=False)
    acid = next(f for f in fragments if f.GetNumAtoms() != len(alkyl_atoms) or f.GetNumAtoms() == mol.GetNumAtoms() - len(alkyl_atoms))
    Chem.SanitizeMol(acid)
    anion = _anion_name(smiles_to_iupac(Chem.MolToSmiles(acid)))

    halogens = halogen_substituents(mol)
    alkyl_name, _ = name_branch(graph, alkyl_carbon, ester_oxygen, halogens, mol=mol)
    free = len(acid.GetSubstructMatches(_FREE_ACID)) - 1
    if free > 0:
        if free not in _HYDROGEN_WORDS:
            raise UnsupportedStructure("too many free acid groups beside the ester")
        return f"{alkyl_name} {_HYDROGEN_WORDS[free]} {anion}"
    return f"{alkyl_name} {anion}"
