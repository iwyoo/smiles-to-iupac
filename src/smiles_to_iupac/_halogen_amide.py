"""Amides of the hypohalous acids, R-NH-X and R2N-X (P-61.3.2.2, P-68.5.3): a halogen atom on nitrogen makes an
amide of an inorganic acid, which outranks the amine and halo classes of the 'N-halogenoamine' name, so the PIN is
'<substituents>hypochlorous amide' with the carbon groups on nitrogen cited without locants."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_mononuclear_prefixes, name_branch

_ACID_WORDS = {9: "hypofluorous", 17: "hypochlorous", 35: "hypobromous", 53: "hypoiodous"}


def _halogen_amide_nitrogen(mol):
    found = [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 7
        and sum(n.GetAtomicNum() in HALOGEN_PREFIXES for n in a.GetNeighbors()) == 1
    ]
    return found[0] if len(found) == 1 else None


def has_halogen_amide_shape(mol):
    return _halogen_amide_nitrogen(mol) is not None


def name_halogen_amide(mol) -> str:
    nitrogen = _halogen_amide_nitrogen(mol)
    if nitrogen is None:
        raise UnsupportedStructure("no nitrogen bearing exactly one halogen atom")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetAtomicNum() not in (6, 7, *HALOGEN_PREFIXES) or a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("only carbon groups and halogens may accompany an amide of a halogen acid")
    if sum(a.GetAtomicNum() == 7 for a in mol.GetAtoms()) != 1 or nitrogen.IsInRing():
        raise UnsupportedStructure("a ring nitrogen or a second nitrogen is not an amide of a halogen acid")
    halogen = next(n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() in HALOGEN_PREFIXES)
    carbons = [n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 6]
    if not carbons or any(
        mol.GetBondBetweenAtoms(nitrogen.GetIdx(), c.GetIdx()).GetBondTypeAsDouble() != 1.0 for c in carbons
    ):
        raise UnsupportedStructure("an amide of a halogen acid needs single-bonded carbon groups on nitrogen")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [
        name_branch(graph, c.GetIdx(), nitrogen.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True)
        for c in carbons
    ]
    return f"{format_mononuclear_prefixes(entries)}{_ACID_WORDS[halogen.GetAtomicNum()]} amide"
