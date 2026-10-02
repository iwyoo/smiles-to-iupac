"""Structural detection (not yet full naming) for furan ortho-fused onto
pyran, both ring oxygens (P-25.3.2.4(c): pyran is base; Blue Book's own
"2H-furo[3,2-b]pyran (PIN)", `tmp/bluebook/P2.txt` ~3998). The fusion
descriptor is verified correct, but this project's own whole-system
numbering for this hexagon+pentagon pair contradicts that "2H" --
`name_furan_pyran_fusion` refuses to guess rather than risk a wrong H."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._pyridine_heterocycle_fusion import _base_bond_letter
from ._two_component_heterocycle_fusion import _local_numbering


def find_furan_pyran_fusion_core(mol):
    """Return (pyran_ring, furan_ring, pyran_o, furan_o, fusion_atoms,
    sp3_carbon) if `mol` is exactly furan ortho-fused to pyran, else
    None."""
    if mol.GetNumAtoms() != 9:
        return None
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None
    sizes = sorted(len(r) for r in atom_rings)
    if sizes != [5, 6]:
        return None
    pyran_ring = next(r for r in atom_rings if len(r) == 6)
    furan_ring = next(r for r in atom_rings if len(r) == 5)

    for idx in furan_ring:
        atom = mol.GetAtomWithIdx(idx)
        if not atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
            return None
    furan_hetero = [idx for idx in furan_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(furan_hetero) != 1 or mol.GetAtomWithIdx(furan_hetero[0]).GetAtomicNum() != 8:
        return None
    furan_o = furan_hetero[0]

    pyran_hetero = [idx for idx in pyran_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(pyran_hetero) != 1 or mol.GetAtomWithIdx(pyran_hetero[0]).GetAtomicNum() != 8:
        return None
    pyran_o = pyran_hetero[0]
    if mol.GetAtomWithIdx(pyran_o).GetTotalNumHs() != 0 or mol.GetAtomWithIdx(pyran_o).GetFormalCharge() != 0:
        return None

    shared = set(pyran_ring) & set(furan_ring)
    if len(shared) != 2 or pyran_o in shared or furan_o in shared:
        return None
    a, b = shared
    if b not in {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()}:
        return None

    sp3 = [
        idx
        for idx in pyran_ring
        if idx not in shared
        and idx != pyran_o
        and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
        and all(bond.GetBondTypeAsDouble() == 1.0 for bond in mol.GetAtomWithIdx(idx).GetBonds())
    ]
    if len(sp3) != 1:
        return None
    for idx in pyran_ring:
        if idx in shared or idx == pyran_o or idx == sp3[0]:
            continue
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0:
            return None
        if not any(bond.GetBondTypeAsDouble() == 2.0 for bond in atom.GetBonds()):
            return None

    if mol.GetNumAtoms() != len(set(pyran_ring) | set(furan_ring)):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    return pyran_ring, furan_ring, pyran_o, furan_o, shared, sp3[0]


def has_furan_pyran_fusion_name(mol) -> bool:
    return find_furan_pyran_fusion_core(mol) is not None


def name_furan_pyran_fusion(mol) -> str:
    core = find_furan_pyran_fusion_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "this two-ring system is not a supported furan + pyran "
            "ortho-fusion (see P-25.3.1.3)"
        )
    pyran_ring, furan_ring, pyran_o, furan_o, fusion_atoms, _sp3_carbon = core
    graph = {atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()}

    base = _base_bond_letter(graph, pyran_ring, pyran_o, fusion_atoms)
    if base is None:
        raise UnsupportedStructure(
            "fusion at pyran's own oxygen is not supported (see P-25.3.1.3)"
        )
    letter, base_numbering = base

    furan_numbering = _local_numbering(graph, furan_ring, furan_o, fusion_atoms)
    if furan_numbering is None:
        raise UnsupportedStructure(
            "the fusion bond does not touch furan's own oxygen-adjacent "
            "carbon (see P-25.3.1.3)"
        )

    base_low, base_high = sorted(fusion_atoms, key=lambda atom: base_numbering[atom])
    citation = f"{furan_numbering[base_low]},{furan_numbering[base_high]}"

    raise UnsupportedStructure(
        f"furo[{citation}-{letter}]pyran's fusion descriptor is known, but "
        "this project's whole-system P-25.3.3.1.1 numbering for a "
        "hexagon+pentagon heterocyclic pair is not yet verified against "
        "the Blue Book's own indicated-hydrogen answer -- "
        "refusing to guess the H locant rather than risk a silently "
        "wrong name"
    )
