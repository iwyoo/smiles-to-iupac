"""Amides of sulfuric and sulfurous acid (P-67.1.2.4.1.1, P-67.1.2.6.1): the acid with one hydroxy group replaced by an
amino group is 'sulfamic acid' (anion 'sulfamate'); with both hydroxy groups replaced it is the 'sulfuric diamide' (or
'sulfurous diamide'), the amino nitrogens being N and N' and each carbon group on them an N prefix.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents
from ._substituents import format_substituent_prefixes, name_branch


def _centre(mol):
    sulfurs = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 16]
    return sulfurs[0] if len(sulfurs) == 1 else None


def _slots(mol, sulfur):
    """(oxo count, [('N', atom) | ('O', atom)] single-bonded ligands) or None."""
    oxo, ligands = 0, []
    for n in sulfur.GetNeighbors():
        order = mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondTypeAsDouble()
        if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and not n.GetFormalCharge() and order == 2.0:
            oxo += 1
        elif n.GetAtomicNum() == 8 and n.GetDegree() == 1 and order == 1.0 and (
            n.GetFormalCharge() == -1 or (not n.GetFormalCharge() and n.GetTotalNumHs() == 1)
        ):
            ligands.append(("O", n))
        elif n.GetAtomicNum() == 7 and order == 1.0 and not n.GetFormalCharge() and not n.IsInRing() and not n.GetIsAromatic():
            ligands.append(("N", n))
        else:
            return None
    return oxo, ligands


def has_sulfuric_amide_shape(mol) -> bool:
    sulfur = _centre(mol)
    if sulfur is None or sulfur.GetFormalCharge() or sulfur.GetIsotope() or len(Chem.GetMolFrags(mol)) != 1:
        return False
    slots = _slots(mol, sulfur)
    if slots is None:
        return False
    oxo, ligands = slots
    kinds = [kind for kind, _ in ligands]
    return oxo in (1, 2) and len(ligands) == 2 and "N" in kinds and sulfur.GetDegree() == oxo + 2


def name_sulfuric_amide(mol) -> str:
    sulfur = _centre(mol)
    oxo, ligands = _slots(mol, sulfur)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    centre = sulfur.GetIdx()
    nitrogens = [n for kind, n in ligands if kind == "N"]
    cited = {}
    for nitrogen in nitrogens:
        entries = []
        for n in graph[nitrogen.GetIdx()]:
            if n == centre:
                continue
            atom = mol.GetAtomWithIdx(n)
            if atom.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n).GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure("only carbon groups on the nitrogen of a sulfuric amide are supported")
            entries.append(name_branch(graph, n, nitrogen.GetIdx(), halogens, mol=mol))
        cited[nitrogen.GetIdx()] = entries
    ordered = sorted(nitrogens, key=lambda a: (-len(cited[a.GetIdx()]), sorted(name for name, _ in cited[a.GetIdx()])))
    locants = ("N", "N'")
    entries = {}
    for locant, nitrogen in zip(locants, ordered):
        for item in cited[nitrogen.GetIdx()]:
            entries.setdefault(locant, []).append(item)
    grouped = group_substituents(entries)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    kinds = [kind for kind, _ in ligands]
    if len(nitrogens) == 2:
        word = "sulfuric diamide" if oxo == 2 else "sulfurous diamide"
        return prefix + word
    if oxo != 2:
        raise UnsupportedStructure("the amide acid of sulfurous acid is not named here")
    anionic = any(n.GetFormalCharge() == -1 for kind, n in ligands if kind == "O")
    return prefix + ("sulfamate" if anionic else "sulfamic acid")
