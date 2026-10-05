"""Chalcogen (S/Se/Te) ring-oxide of a bare fused ring system, named by P-62.5's functional-class "oxide" with the
locant of the chalcogen in the numbering of the fused parent (P-25.3.3)."""

from rdkit import Chem

from ._fusion_name import Context, fused_ring_system_name, fusion_name, system_numbering_options
from ._fused_numbering import _locant_key

_CHALCOGENS = {16, 34, 52}


def _reduced_ring(mol, oxide_oxygen_idx):
    reduced = Chem.RWMol(mol)
    reduced.RemoveAtom(oxide_oxygen_idx)
    reduced_mol = reduced.GetMol()
    Chem.SanitizeMol(reduced_mol)
    return reduced_mol


def _find_fused_ring_oxide(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return None
    ring_atoms = {idx for ring in ring_info.AtomRings() for idx in ring}
    candidates = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetIdx() in ring_atoms
        and atom.GetAtomicNum() in _CHALCOGENS
        and atom.GetFormalCharge() == 0
        and atom.GetIsotope() == 0
    ]
    if len(candidates) != 1:
        return None
    (heteroatom,) = candidates
    oxide_neighbors = [
        n
        for n in heteroatom.GetNeighbors()
        if n.GetIdx() not in ring_atoms
        and n.GetAtomicNum() == 8
        and n.GetFormalCharge() == 0
        and n.GetIsotope() == 0
        and n.GetDegree() == 1
    ]
    if len(oxide_neighbors) != 1:
        return None
    (oxide_oxygen,) = oxide_neighbors
    if mol.GetBondBetweenAtoms(heteroatom.GetIdx(), oxide_oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
        return None
    return heteroatom.GetIdx(), oxide_oxygen.GetIdx()


def _match(mol):
    found = _find_fused_ring_oxide(mol)
    if found is None:
        return None
    heteroatom, oxide_oxygen = found
    reduced = _reduced_ring(mol, oxide_oxygen)
    base_name = fused_ring_system_name(reduced)
    if base_name is None:
        return None
    ctx = Context(reduced)
    name, root = fusion_name(reduced)
    options = system_numbering_options(ctx, name, root)
    locant = min((n[heteroatom - (1 if heteroatom > oxide_oxygen else 0)] for n in options), key=_locant_key)
    return base_name, locant


def has_fused_hetero_ring_oxide_shape(mol) -> bool:
    return _match(mol) is not None


def name_fused_hetero_ring_oxide(mol) -> str:
    base_name, locant = _match(mol)
    return f"{base_name} {locant}-oxide"
