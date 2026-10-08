"""Halogen oxo groups -XO, -XO2, -XO3 as compulsory prefixes (P-61.3.2.3, P-67.1.4.5): 'chlorosyl', 'chloryl' and
'perchloryl' (and the bromo, iodo and fluoro analogues) on a saturated acyclic chain or a plain benzene ring. Like
'azido' they have no suffix form, so the parent is chosen by the ordinary hydrocarbon rules (P-44.1.2.2)."""

from rdkit import Chem

from ._acyclic import name_from_carbon_graph
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    halogen_substituents,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
)
from ._hetero_prefixes import halogen_oxo_prefix
from ._multiplicative_text import enclose
from ._substituents import name_branch


def _oxo_halogens(mol):
    """{halogen idx: prefix} for every halogen carrying terminal oxygens and one carbon."""
    found = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in HALOGEN_PREFIXES or atom.GetDegree() < 2:
            continue
        carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbons) != 1:
            continue
        prefix = halogen_oxo_prefix(mol, atom.GetIdx(), carbons[0].GetIdx())
        if prefix is not None:
            found[atom.GetIdx()] = prefix if prefix.startswith(("chlor", "brom", "iod", "fluor", "per")) else enclose(prefix)
    return found


def has_halogen_oxo_shape(mol) -> bool:
    return bool(_oxo_halogens(mol))


def name_halogen_oxo(mol) -> str:
    groups = _oxo_halogens(mol)
    if not groups:
        raise UnsupportedStructure("no halogen oxo group found")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    owned = set(groups)
    for idx in groups:
        owned.update(n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetAtomicNum() in (8, 16, 34, 52))
    ring_atoms = frozenset()
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1 and is_plain_benzene_ring(mol, set(ring_info.AtomRings()[0])):
        ring_atoms = frozenset(ring_info.AtomRings()[0])
    elif ring_info.NumRings():
        raise UnsupportedStructure("rings other than one plain benzene ring are not supported here")
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in owned or idx in ring_atoms:
            continue
        if atom.GetIsotope() or atom.GetFormalCharge():
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() != 6 and not (atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetDegree() == 1):
            raise UnsupportedStructure("only carbon and halogen substituents may accompany a halogen oxo group")
        if atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic atoms outside one plain benzene ring are not supported yet")
    if any(a not in owned and b not in owned and a not in ring_atoms for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure("unsaturation is not supported by this module")

    terminals = {**halogen_substituents(mol), **groups}
    graph = adjacency(mol)
    if ring_atoms:
        attachment = ring_chain_attachment(graph, set(ring_atoms), set())
        if attachment is None:
            raise UnsupportedStructure("a benzene ring with several substituents is not supported here")
        ring_atom, root = attachment
        name, compound = name_branch(graph, root, ring_atom, terminals, mol=mol)
        return f"{enclose(name) if compound else name}benzene"
    return name_from_carbon_graph(graph, carbon_adjacency(mol), terminals, mol=mol)
