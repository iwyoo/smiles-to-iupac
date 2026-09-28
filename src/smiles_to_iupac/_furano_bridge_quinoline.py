"""P-25.4.2.1.5 heterocyclic-bridge citation for an intact furan ring
whose own adjacent C2/C3 atoms each form a new bond to a non-adjacent
pair of `_pyridine_bicyclic_fusion.py` base atoms (benzo[g]quinoline's
meso locants), rather than sharing atoms with the base. Verified against
the Blue Book's own diagram (`tmp/bluebook/P2.pdf` p.106, cached at
`tmp/bluebook/p25_bridge_examples/page106-106.png`): `10,5-[2,3]furano-
benzo[g]quinoline` -- no PubChem-registered structure exists for this
worked example, so this is deductively verified against the primary
source only, not cross-checked against a registry entry.
"""

from rdkit import Chem

from ._parent_hydride_stripping import strip_substituents
from ._pyridine_bicyclic_fusion import has_pyridine_bicyclic_fusion_name, name_pyridine_bicyclic_fusion
from ._quinoline_bicyclic_numbering import peripheral_numbering


def _find_furan_bridge(mol):
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetIsAromatic() for a in atoms):
            continue
        o_atoms = [a.GetIdx() for a in atoms if a.GetAtomicNum() == 8]
        c_atoms = [a.GetIdx() for a in atoms if a.GetAtomicNum() == 6]
        if len(o_atoms) != 1 or len(c_atoms) != 4:
            continue
        ring_set = set(ring)
        external = [
            (a.GetIdx(), n.GetIdx())
            for a in atoms
            for n in a.GetNeighbors()
            if n.GetIdx() not in ring_set
        ]
        if len(external) != 2:
            continue
        (f1, b1), (f2, b2) = external
        if mol.GetBondBetweenAtoms(f1, f2) is None:
            continue
        o_idx = o_atoms[0]
        f1_adj_o = mol.GetBondBetweenAtoms(f1, o_idx) is not None
        f2_adj_o = mol.GetBondBetweenAtoms(f2, o_idx) is not None
        if f1_adj_o == f2_adj_o:
            continue
        if f1_adj_o:
            return ring_set, b1, b2
        return ring_set, b2, b1
    return None


def find_furano_bridge_quinoline(mol):
    """Return (base_name, loc_c2, loc_c3) or None."""
    bridge = _find_furan_bridge(mol)
    if bridge is None:
        return None
    furan_ring, base_c2, base_c3 = bridge

    # A genuine ortho-fused (not bridged) 5-ring elsewhere in the molecule
    # can also match `_find_furan_bridge`'s shape; stripping its shared
    # atoms then breaks a neighboring ring, which surfaces as any of
    # several RDKit sanitize exception subclasses -- all mean "not this".
    try:
        stripped, old_to_new = strip_substituents(mol, furan_ring)
    except Chem.rdchem.MolSanitizeException:
        return None
    if stripped is None:
        return None
    if not has_pyridine_bicyclic_fusion_name(stripped):
        return None

    numbering = peripheral_numbering(stripped)
    if numbering is None:
        return None
    loc_c2 = numbering.get(old_to_new.get(base_c2))
    loc_c3 = numbering.get(old_to_new.get(base_c3))
    if loc_c2 is None or loc_c3 is None or not loc_c2.isdigit() or not loc_c3.isdigit():
        return None

    return name_pyridine_bicyclic_fusion(stripped), loc_c2, loc_c3


def has_furano_bridge_quinoline_name(mol) -> bool:
    return find_furano_bridge_quinoline(mol) is not None


def name_furano_bridge_quinoline(mol) -> str:
    base_name, loc_c2, loc_c3 = find_furano_bridge_quinoline(mol)
    return f"{loc_c2},{loc_c3}-[2,3]furano{base_name}"
