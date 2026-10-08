"""Anions of phosphorus, boron and sulfur oxoacids and their esters (P-72.2.2.2.1.1,
P-72.2.2.2.1.2): phosphinate/phosphonate/phosphate and their 'ite' analogues take
the acid ending 'ate'/'ite'; acid esters of inorganic acids are named as 'hydrogen
salts' with the ester groups cited first, then 'hydrogen' (P-65.6.2.3, P-65.6.3.3.5).
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._numerals import multiplying_prefix
from ._substituents import format_mononuclear_prefixes, name_branch

_PHOSPHORUS = {
    (5, 0): "phosphate",
    (5, 1): "phosphonate",
    (5, 2): "phosphinate",
    (3, 0): "phosphite",
    (3, 1): "phosphonite",
    (3, 2): "phosphinite",
}
_BORON = {0: "borate", 1: "boronate", 2: "borinate"}


_HEAVY_STEMS = {33: "ars", 51: "stib"}


def _pnictogen_anion(z, key):
    """'phosphonate', 'arsonate', 'stiborate' ...: the phosphorus anion name with the stem of the element."""
    name = _PHOSPHORUS[key]
    if z == 15:
        return name
    stem = _HEAVY_STEMS[z]
    return stem + ("or" + name[len("phosph"):] if key[1] == 0 else name[len("phosph"):])


def oxoacid_center(mol):
    """The single P or oxygen-only S atom bearing the anionic oxygen, else None."""
    anions = [a for a in mol.GetAtoms() if a.GetFormalCharge() < 0]
    if not anions or any(a.GetAtomicNum() != 8 or a.GetDegree() != 1 or a.GetFormalCharge() != -1 for a in anions):
        return None
    hosts = {a.GetNeighbors()[0].GetIdx() for a in anions}
    if len(hosts) != 1:
        return None
    center = mol.GetAtomWithIdx(next(iter(hosts)))
    if center.GetAtomicNum() in (5, 15, 33, 51):
        return center
    if center.GetAtomicNum() == 16 and not any(n.GetAtomicNum() == 6 for n in center.GetNeighbors()):
        return center
    return None


def _pseudohalide_infix_anion(mol, center):
    """'hydrogen borocyanatidate' (P-67.1.3.1): an -OC#N group on the acid centre is the infix 'cyanatid', so the anion
    takes the name of the acid with a replaced hydroxy group."""
    cyanates = [
        n
        for n in center.GetNeighbors()
        if n.GetAtomicNum() == 8
        and n.GetDegree() == 2
        and any(
            m.GetAtomicNum() == 6
            and m.GetDegree() == 2
            and any(x.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(m.GetIdx(), x.GetIdx()).GetBondTypeAsDouble() == 3.0 for x in m.GetNeighbors())
            for m in n.GetNeighbors()
            if m.GetIdx() != center.GetIdx()
        )
    ]
    if not cyanates or len(mol.GetAtoms()) != len(center.GetNeighbors()) + 1 + 2 * len(cyanates):
        return None
    from ._acid_derivatives import anion_name
    from .core import smiles_to_iupac

    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        if atom.GetFormalCharge() == -1:
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
    neutral = editable.GetMol()
    Chem.SanitizeMol(neutral)
    acid = smiles_to_iupac(Chem.MolToSmiles(neutral))
    if not acid.endswith(" acid"):
        return None
    hydroxyls = sum(
        1
        for n in center.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1 and not n.GetFormalCharge()
    )
    base = anion_name(acid)
    if not hydroxyls:
        return base
    return ("hydrogen" if hydroxyls == 1 else multiplying_prefix(hydroxyls) + "hydrogen") + " " + base


def name_oxoacid_anion(mol, center):
    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("a multi-fragment structure is not named by one oxoacid anion")
    if any(s.specified == Chem.StereoSpecified.Specified for s in Chem.FindPotentialStereo(mol)):
        raise UnsupportedStructure("the stereochemistry of an oxoacid ester anion is not cited yet")
    if len(set(a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() > 0)) or len(mol.GetAtoms()) < 2:
        raise UnsupportedStructure("cationic atoms beside an oxoacid anion are not supported here")
    infix = _pseudohalide_infix_anion(mol, center)
    if infix is not None:
        return infix
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    carbons, oxo, hydroxyls, esters = [], 0, 0, []
    for n in center.GetNeighbors():
        z = n.GetAtomicNum()
        bond = mol.GetBondBetweenAtoms(center.GetIdx(), n.GetIdx()).GetBondTypeAsDouble()
        if z == 6 and bond == 1.0:
            carbons.append(n.GetIdx())
        elif z == 8 and bond == 2.0 and n.GetDegree() == 1:
            oxo += 1
        elif z == 8 and n.GetFormalCharge() == -1:
            continue
        elif z == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1:
            hydroxyls += 1
        elif z == 8 and n.GetDegree() == 2:
            other = next(x for x in n.GetNeighbors() if x.GetIdx() != center.GetIdx())
            if other.GetAtomicNum() != 6:
                raise UnsupportedStructure("this oxoacid ester oxygen is not bonded to carbon")
            esters.append(name_branch(graph, other.GetIdx(), n.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True))
        else:
            raise UnsupportedStructure("this oxoacid substituent is not supported yet")
    if any(a.GetIdx() != center.GetIdx() and a.GetAtomicNum() in (5, 15, 16) for a in mol.GetAtoms()):
        raise UnsupportedStructure("several phosphorus atoms are not supported here")
    if center.GetAtomicNum() == 5:
        if oxo or len(carbons) > 2:
            raise UnsupportedStructure("this boron acid pattern is not supported yet")
        entries = [
            name_branch(graph, c, center.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True) for c in carbons
        ]
        parent = (format_mononuclear_prefixes(entries) if entries else "") + _BORON[len(carbons)]
    elif center.GetAtomicNum() in (15, 33, 51):
        key = (5 if oxo else 3, len(carbons))
        if key not in _PHOSPHORUS:
            raise UnsupportedStructure("this phosphorus oxoacid pattern is not supported yet")
        entries = [
            name_branch(graph, c, center.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True) for c in carbons
        ]
        parent = (format_mononuclear_prefixes(entries) if entries else "") + _pnictogen_anion(center.GetAtomicNum(), key)
    else:
        if oxo == 2:
            parent = "sulfate"
        elif oxo == 1:
            parent = "sulfite"
        else:
            raise UnsupportedStructure("this sulfur oxoacid pattern is not supported yet")
    words = []
    if esters:
        words.append(_ester_words(esters))
    if hydroxyls:
        words.append(("hydrogen" if hydroxyls == 1 else multiplying_prefix(hydroxyls) + "hydrogen"))
    words.append(parent)
    return " ".join(words)


def _ester_words(esters):
    groups = {}
    for name, compound in esters:
        groups.setdefault((name, compound), 0)
        groups[(name, compound)] += 1
    parts = []
    for (name, compound), count in sorted(groups.items(), key=lambda item: item[0][0]):
        text = f"({name})" if compound else name
        parts.append(text if count == 1 else multiplying_prefix(count, compound=compound) + text)
    return " ".join(parts)
