"""Didehydro derivatives of mancude heteromonocycles whose numbering is fixed by the hetero atoms (P-31.2.2, P-44.4.1.10.2):
the ring triple bond of pyridyne is expressed by the prefix 'didehydro' on the parent hydride, 3,4-didehydropyridine; the
lowest locants of the pair go to the bond (P-31.1.4.2.4)."""

from rdkit import Chem

from ._common import adjacency
from ._ring_assembly_chain import _NON_NH_ROLE_SEQUENCES, _hetero_ring_alignments, _match_hetero_ring_parent


def _parent(mol):
    """(parent mol, ring atoms, triple bond atoms) when `mol` is one unsubstituted ring with one ring C#C, else None."""
    ring_info = mol.GetRingInfo()
    if len(Chem.GetMolFrags(mol)) != 1 or ring_info.NumRings() != 1 or len(ring_info.AtomRings()[0]) != mol.GetNumAtoms():
        return None
    triples = [b for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 3.0]
    if len(triples) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    editable = Chem.RWMol(mol)
    a, b = triples[0].GetBeginAtomIdx(), triples[0].GetEndAtomIdx()
    editable.GetBondBetweenAtoms(a, b).SetBondType(Chem.BondType.DOUBLE)
    parent = editable.GetMol()
    try:
        Chem.SanitizeMol(parent)
    except Exception:
        return None
    if not all(x.GetIsAromatic() for x in parent.GetAtoms()):
        return None
    return parent, list(ring_info.AtomRings()[0]), (a, b)


def has_hetaryne_shape(mol) -> bool:
    found = _parent(mol)
    if found is None:
        return False
    parent, ring, _ = found
    name = _match_hetero_ring_parent(parent, adjacency(parent), ring)
    return name in _NON_NH_ROLE_SEQUENCES


def name_hetaryne(mol) -> str:
    parent, ring, (a, b) = _parent(mol)
    graph = adjacency(parent)
    name = _match_hetero_ring_parent(parent, graph, ring)
    best = min(
        tuple(sorted((alignment[a], alignment[b]))) for alignment in _hetero_ring_alignments(parent, graph, ring, name)
    )
    return f"{best[0]},{best[1]}-didehydro{name}"
