"""P-25.4.2.1.5 heterocyclic-bridge citation for an intact pyran ring
whose own para C2/C5 atoms (C2 adjacent to O) each bond to a
`_pyridine_bicyclic_fusion.py` base; 'epi' distinguishes the bridge
prefix from pyran's identical fusion prefix. Verified deductively
against `tmp/bluebook/p25_bridge_examples/page112-112.png` (no PubChem
structure exists for this worked example)."""

from rdkit import Chem

from ._parent_hydride_stripping import strip_substituents
from ._pyridine_bicyclic_fusion import has_pyridine_bicyclic_fusion_name, name_pyridine_bicyclic_fusion
from ._quinoline_bicyclic_numbering import peripheral_numbering


def _locant_sort_key(locant):
    return (int(locant.rstrip("a")), 1 if locant.endswith("a") else 0)


def _find_pyran_bridge(mol):
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        o_atoms = [a.GetIdx() for a in atoms if a.GetAtomicNum() == 8]
        c_atoms = [a.GetIdx() for a in atoms if a.GetAtomicNum() == 6]
        if len(o_atoms) != 1 or len(c_atoms) != 5:
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
        if mol.GetBondBetweenAtoms(f1, f2) is not None:
            continue
        o_idx = o_atoms[0]
        f1_adj_o = mol.GetBondBetweenAtoms(f1, o_idx) is not None
        f2_adj_o = mol.GetBondBetweenAtoms(f2, o_idx) is not None
        if f1_adj_o == f2_adj_o:
            continue
        sp3 = [
            a.GetIdx()
            for a in atoms
            if a.GetIdx() not in (f1, f2)
            and a.GetIdx() != o_idx
            and all(b.GetBondType() != Chem.BondType.DOUBLE for b in a.GetBonds())
        ]
        if len(sp3) != 1:
            continue
        if f1_adj_o:
            return ring_set, o_idx, f1, b1, f2, b2, sp3[0]
        return ring_set, o_idx, f2, b2, f1, b1, sp3[0]
    return None


def find_pyrano_bridge_quinoline(mol):
    """Return (base_name, loc_c2, loc_c5, indicated_h_locant) or None."""
    bridge = _find_pyran_bridge(mol)
    if bridge is None:
        return None
    pyran_ring, o_idx, c2, base_c2, c5, base_c5, sp3_idx = bridge

    try:
        stripped, old_to_new = strip_substituents(mol, pyran_ring)
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
    loc_c5 = numbering.get(old_to_new.get(base_c5))
    if loc_c2 is None or loc_c5 is None or not loc_c2.isdigit() or not loc_c5.isdigit():
        return None

    highest = max(numbering.values(), key=_locant_sort_key)
    bridge_start = int(highest.rstrip("a")) + 1

    # P-25.4.4: bridge numbering starts at the atom bonded to the
    # bridgehead with the higher locant, then walks the whole ring in
    # whichever direction gives the heteroatom the lower locant.
    start = c2 if int(loc_c2) > int(loc_c5) else c5
    neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(start).GetNeighbors() if n.GetIdx() in pyran_ring]

    def _walk(first_step):
        order = [start]
        prev, cur = start, first_step
        while len(order) < len(pyran_ring):
            order.append(cur)
            nxt = next(
                n.GetIdx()
                for n in mol.GetAtomWithIdx(cur).GetNeighbors()
                if n.GetIdx() in pyran_ring and n.GetIdx() != prev
            )
            prev, cur = cur, nxt
        return order

    candidates = [_walk(n) for n in neighbors]
    order = min(candidates, key=lambda seq: seq.index(o_idx))
    bridge_locants = {atom: str(bridge_start + i) for i, atom in enumerate(order)}
    indicated_h = bridge_locants[sp3_idx]

    return name_pyridine_bicyclic_fusion(stripped), loc_c2, loc_c5, indicated_h


def has_pyrano_bridge_quinoline_name(mol) -> bool:
    return find_pyrano_bridge_quinoline(mol) is not None


def name_pyrano_bridge_quinoline(mol) -> str:
    base_name, loc_c2, loc_c5, indicated_h = find_pyrano_bridge_quinoline(mol)
    return f"{indicated_h}H-{loc_c2},{loc_c5}-[2,5]epipyrano{base_name}"
