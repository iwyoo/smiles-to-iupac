"""Esters of the halogen oxoacids, R-E-X and R-O-XO(n) (P-61.3.2.2, P-67.1.3.2): '<R> hypochlorite', 'chlorite',
'chlorate', 'perchlorate' (and the fluorine, bromine and iodine analogues), with 'thio', 'seleno' or 'telluro'
inserted when the ester atom is a heavier chalcogen: 'phenyl hypofluorite', 'methyl thiohypochlorite',
'3-oxobutyl bromate'. A halogen on a heteroatom is an acid ester, which outranks the substitutive halo name.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._hetero_prefixes import halogen_oxo_prefix
from ._phosphanone import _SENIOR
from ._substituents import name_branch

_STEMS = {9: "fluor", 17: "chlor", 35: "brom", 53: "iod"}
_INFIXES = {8: "", 16: "thio", 34: "seleno", 52: "telluro"}


def _ester_halogens(mol):
    """[(halogen, ester atom, carbon, number of oxo oxygens)] for every R-E-X(O)n group."""
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in HALOGEN_PREFIXES:
            continue
        hosts = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in _INFIXES and n.GetDegree() == 2]
        if len(hosts) != 1:
            continue
        host = hosts[0]
        carbons = [n for n in host.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbons) != 1:
            continue
        oxo = [n for n in atom.GetNeighbors() if n.GetIdx() != host.GetIdx()]
        if len(oxo) > 3 or any(n.GetAtomicNum() != 8 or n.GetDegree() != 1 for n in oxo):
            continue
        if oxo and halogen_oxo_prefix(mol, atom.GetIdx(), host.GetIdx()) is None:
            continue
        if not oxo and atom.GetFormalCharge():
            continue
        found.append((atom, host, carbons[0], len(oxo)))
    return found


def has_halogen_acid_ester_shape(mol) -> bool:
    return len(_ester_halogens(mol)) == 1


def name_halogen_acid_ester(mol) -> str:
    groups = _ester_halogens(mol)
    if len(groups) != 1:
        raise UnsupportedStructure("exactly one halogen oxoacid ester group is supported")
    halogen, host, carbon, oxygens = groups[0]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if oxygens and host.GetAtomicNum() != 8:
        raise UnsupportedStructure("a chalcogen analogue of a halogen oxoacid with oxygens is not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR):
        raise UnsupportedStructure("an acid, ester, amide, nitrile or aldehyde group outranks the halogen acid ester")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    name, _ = name_branch(graph, carbon.GetIdx(), host.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
    stem = _STEMS[halogen.GetAtomicNum()]
    acid = {0: f"hypo{stem}ite", 1: f"{stem}ite", 2: f"{stem}ate", 3: f"per{stem}ate"}[oxygens]
    return f"{name} {_INFIXES[host.GetAtomicNum()]}{acid}"
