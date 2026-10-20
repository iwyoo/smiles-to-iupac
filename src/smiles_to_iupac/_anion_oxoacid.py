"""Anions of phosphorus, boron and sulfur oxoacids and their esters (P-72.2.2.2.1.1,
P-72.2.2.2.1.2): phosphinate/phosphonate/phosphate and their 'ite' analogues take
the acid ending 'ate'/'ite'; acid esters of inorganic acids are named as 'hydrogen
salts' with the ester groups cited first, then 'hydrogen' (P-65.6.2.3, P-65.6.3.3.5).
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._hetero_prefixes import _ONIUM_PREFIX_STEMS
from ._numerals import multiplying_prefix
from ._phosphate import citing_stereo, format_ester_entries
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


def _alkoxide_host(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    return atom.GetAtomicNum() == 6 and not any(
        b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16, 34, 52) for b in atom.GetBonds()
    )


def oxoacid_center(mol):
    """The single P or oxygen-only S atom bearing the anionic oxygen, else None. Alkoxide groups elsewhere are cited as
    'oxido' prefixes: the phosphorus acid anion is the senior parent (P-72.7(e))."""
    anions = [a for a in mol.GetAtoms() if a.GetFormalCharge() < 0]
    if not anions or any(a.GetAtomicNum() != 8 or a.GetDegree() != 1 or a.GetFormalCharge() != -1 for a in anions):
        return None
    hosts = {a.GetNeighbors()[0].GetIdx() for a in anions}
    if len(hosts) != 1:
        centers = {h for h in hosts if mol.GetAtomWithIdx(h).GetAtomicNum() in (15, 33, 51)}
        if len(centers) != 1 or not all(_alkoxide_host(mol, h) for h in hosts - centers):
            return None
        hosts = centers
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
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge() > 0]
    anions = [a for a in mol.GetAtoms() if a.GetFormalCharge() < 0]
    ammonium_only = all(
        (a.GetAtomicNum() == 7 or a.GetAtomicNum() in _ONIUM_PREFIX_STEMS) and a.GetFormalCharge() == 1 and not a.GetIsotope()
        for a in cations
    )
    if len(mol.GetAtoms()) < 2 or (cations and not (ammonium_only and len(cations) == len(anions))):
        raise UnsupportedStructure("cationic atoms beside an oxoacid anion are not supported here")
    infix = _pseudohalide_infix_anion(mol, center)
    if infix is not None:
        return infix
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    with citing_stereo(mol):
        carbons, oxo, hydroxyls, esters, anionic = [], 0, 0, [], 0
        for n in center.GetNeighbors():
            z = n.GetAtomicNum()
            bond = mol.GetBondBetweenAtoms(center.GetIdx(), n.GetIdx()).GetBondTypeAsDouble()
            if z == 6 and bond == 1.0:
                carbons.append(n.GetIdx())
            elif z == 8 and bond == 2.0 and n.GetDegree() == 1:
                oxo += 1
            elif z == 8 and n.GetFormalCharge() == -1:
                anionic += 1
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
        if center.GetAtomicNum() == 5:
            if oxo or len(carbons) > 2:
                raise UnsupportedStructure("this boron acid pattern is not supported yet")
            entries = [
                name_branch(graph, c, center.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True) for c in carbons
            ]
            parent = (format_mononuclear_prefixes(entries) if entries else "") + _BORON[len(carbons)]
        elif center.GetAtomicNum() in (15, 33, 51):
            key = (5 if oxo else 3, 3 - anionic - hydroxyls - len(esters))
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


def multiplied_phosphonate_name(mol):
    """'methylenebis(hydrogen phosphonate)': identical phosphonate groups, each with the same number of anionic oxygens,
    joined by a linker take the multiplicative name of the neutral acid with 'phosphonic acid' turned into the anion
    (P-72.2.2.2, P-67.1.3.1)."""
    import re

    from .core import smiles_to_iupac

    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() > 0 for a in mol.GetAtoms()):
        return None
    anions = [a for a in mol.GetAtoms() if a.GetFormalCharge() < 0]
    if not anions or any(a.GetAtomicNum() != 8 or a.GetDegree() != 1 or a.GetFormalCharge() != -1 for a in anions):
        return None
    hosts = {}
    for a in anions:
        hosts.setdefault(a.GetNeighbors()[0].GetIdx(), []).append(a.GetIdx())
    phosphorus = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    if len(hosts) < 2 or len(phosphorus) != len(hosts) or len({len(o) for o in hosts.values()}) != 1:
        return None
    if any(mol.GetAtomWithIdx(h).GetAtomicNum() != 15 for h in hosts):
        return None
    if any(len(h.GetNeighbors()) != 4 for h in phosphorus) or any(
        sum(1 for n in mol.GetAtomWithIdx(h).GetNeighbors() if n.GetAtomicNum() == 6) != 1 for h in hosts
    ):
        return None
    editable = Chem.RWMol(mol)
    for a in anions:
        atom = editable.GetAtomWithIdx(a.GetIdx())
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    neutral = editable.GetMol()
    Chem.SanitizeMol(neutral)
    acid = smiles_to_iupac(Chem.MolToSmiles(neutral))
    found = re.fullmatch(r"(?P<head>.+)(?P<multiplier>bis|tris|tetrakis)\(phosphonic acid\)", acid)
    if found is None:
        return None
    unit = "phosphonate" if len(next(iter(hosts.values()))) == 2 else "hydrogen phosphonate"
    return f"{found.group('head')}{found.group('multiplier')}({unit})"


def _ester_words(esters):
    return format_ester_entries(esters)
