"""Hydroxylamine as a functional parent carrying an acid or amide suffix on the oxygen (P-68.3.1.1.1.4):
H2N-O-SO2-OH is 'hydroxylamine-O-sulfonic acid', H2N-O-CO-OH 'hydroxylamine-O-carboxylic acid', H2N-O-CO-NH2
'hydroxylamine-O-carboxamide'; alkyl groups on the nitrogen are cited as N-prefixes."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_substituent_prefixes, name_branch

_SUFFIXES = (
    ("[CX3](=O)[OX2H1]", "carboxylic acid"),
    ("[CX3](=O)[NX3H2]", "carboxamide"),
    ("[SX4](=O)(=O)[OX2H1]", "sulfonic acid"),
    ("[SX4](=O)(=O)[NX3H2]", "sulfonamide"),
    ("[SX3](=O)[OX2H1]", "sulfinic acid"),
)


def _match(mol):
    """(nitrogen, oxygen, acyl atom, suffix) of an H2N-O-acyl or RHN-O-acyl group with nothing else in the molecule."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    for nitrogen in (a for a in mol.GetAtoms() if a.GetAtomicNum() == 7):
        oxygens = [n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 2]
        if len(oxygens) != 1 or nitrogen.IsInRing():
            continue
        oxygen = oxygens[0]
        acyl = next(n for n in oxygen.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx())
        others = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != oxygen.GetIdx()]
        if any(n.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in others):
            continue
        for smarts, suffix in _SUFFIXES:
            for hit in mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
                group = set(hit)
                if hit[0] != acyl.GetIdx() or not {n.GetIdx() for n in acyl.GetNeighbors()} <= group | {oxygen.GetIdx()}:
                    continue
                rest = {a.GetIdx() for a in mol.GetAtoms()} - group - {nitrogen.GetIdx(), oxygen.GetIdx()}
                if all(mol.GetAtomWithIdx(i).GetAtomicNum() in (6, 9, 17, 35, 53) for i in rest):
                    return nitrogen, oxygen, acyl, suffix
    return None


def has_hydroxylamine_acid_shape(mol) -> bool:
    return _match(mol) is not None


def name_hydroxylamine_acid(mol) -> str:
    found = _match(mol)
    if found is None:
        raise UnsupportedStructure("not a hydroxylamine-O-acid")
    nitrogen, oxygen, _, suffix = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    grouped = {}
    for n in nitrogen.GetNeighbors():
        if n.GetIdx() == oxygen.GetIdx():
            continue
        name, compound = name_branch(graph, n.GetIdx(), nitrogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append("N")
    return f"{format_substituent_prefixes(grouped) if grouped else ''}hydroxylamine-O-{suffix}"
