"""Acids and acid halides of the pentavalent nitrogen acids substituted by organyl groups (P-67.1.1.2, P-67.1.2.5.1):
R-N(O)(OH)2 'methylazonic acid', R2N(O)(OH) 'dimethylazinic acid', R-N(O)X2 'methylazonic dichloride' and
R2N(O)X 'dimethylazinic chloride'. The N(=O) group is the zwitterion N(+)-O(-) in the structure RDKit builds."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._hetero_prefixes import pentavalent_nitrogen_oxo
from ._phosphonic_acid import acid_prefixes, senior_acid_present
from ._substituents import format_mononuclear_prefixes, name_branch

_HALIDE = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_MULTIPLIER = {1: "", 2: "di"}


def _find(mol):
    """(nitrogen, organyl carbons, ligands, 'acid' | 'halide') for a molecule built on one pentavalent nitrogen with
    one oxo group, one or two organyl groups and hydroxy groups or identical halogens for the rest."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    centres = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7 and a.GetDegree() == 4 and not a.IsInRing()]
    if len(centres) != 1:
        return None
    nitrogen = centres[0]
    oxo = [n for n in nitrogen.GetNeighbors() if pentavalent_nitrogen_oxo(mol, nitrogen.GetIdx(), n.GetIdx())]
    if len(oxo) != 1:
        return None
    rest = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != oxo[0].GetIdx()]
    if any(mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in rest):
        return None
    carbons = [n for n in rest if n.GetAtomicNum() == 6]
    ligands = [n for n in rest if n.GetAtomicNum() != 6]
    if len(carbons) not in (1, 2) or not ligands:
        return None
    if all(n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1 and not n.GetFormalCharge() for n in ligands):
        return nitrogen, carbons, ligands, "acid"
    elements = {n.GetAtomicNum() for n in ligands}
    if len(elements) == 1 and elements <= set(HALOGEN_PREFIXES) and all(n.GetDegree() == 1 and not n.GetFormalCharge() for n in ligands):
        return nitrogen, carbons, ligands, "halide"
    return None


def has_azonic_acid_shape(mol) -> bool:
    return _find(mol) is not None


def name_azonic_acid(mol) -> str:
    found = _find(mol)
    if found is None:
        raise UnsupportedStructure("no azonic or azinic acid shape")
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified structures are not supported yet")
    if senior_acid_present(mol):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the azonic acid")
    nitrogen, carbons, ligands, kind = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    with acid_prefixes():
        entries = [name_branch(graph, c.GetIdx(), nitrogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True) for c in carbons]
    organyl = format_mononuclear_prefixes(entries)
    acid = "azonic" if len(carbons) == 1 else "azinic"
    if kind == "acid":
        return f"{organyl}{acid} acid"
    return f"{organyl}{acid} {_MULTIPLIER[len(ligands)]}{_HALIDE[ligands[0].GetAtomicNum()]}"
