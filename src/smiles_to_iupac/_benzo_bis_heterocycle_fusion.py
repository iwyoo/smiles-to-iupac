"""Multiparent fusion naming (P-25.3.4.1.3) for a three-ring system: one
plain benzo ring bridging two *identical* five-membered heteromonocycles
(furan or thiophene), each ortho-fused to the benzo ring at a separate
bond, per the IUPAC 2013 Recommendations ("the Blue Book").

This is a multiparent name, not the simpler "attached component(s) on one
parent" shape it might look like at a glance: P-25.3.2.4(a) makes each
furan/thiophene (containing O/S) senior to the all-carbon benzo ring, and
since there are *two* of them, P-25.3.4.1.2's multiplying-prefix rule
applies to the *parent* ("di" + the parent name, e.g. "dithiophene") while
the single benzo ring is the *attached*, bridging component - directly
analogous to the Blue Book's own P-25.3.4.1.3 worked example, `cyclopenta
[1,2-b:1,5-b']difuran` (attached = cyclopenta, two furan parents).
Confirmed against real registered examples below rather than assumed.

Two numbering problems, reusing existing machinery where possible:

- Each parent's own fusion-bond letter: identical to
  `_pyridine_heterocycle_fusion.py`'s `_base_bond_letter` (a ring with one
  fixed-locant-1 heteroatom, direction chosen to give the fusion bond the
  earliest letter) - reused directly, since that function's logic never
  actually depended on the fixed atom being nitrogen specifically.
- The benzo ring's own numbering must satisfy *two* fusion bonds at once
  (one per parent connection), unlike every prior 2-component module here
  which only ever needed one. `_benzo_numbering` tries all 12 (start,
  direction) combinations (benzo has no heteroatom to fix locant 1, so
  every rotation and both directions are candidates), and for each,
  determines which of the two fusion bonds gets the lower locant set (that
  one is cited first, unprimed; the other is primed, per P-25.3.4.1.1's
  "second [parent occurrence] gets primed locants") - then picks whichever
  candidate makes the *first-cited* (unprimed) locant set as low as
  possible, the natural extension of P-25.3.1.3's own "as low as
  possible" rule to two simultaneous fusion bonds.

Real registered examples confirming this family and scope: benzo[1,2-b:
4,5-b']dithiophene (PubChem CID 11106168, a well-known organic-
semiconductor scaffold, "BDT") and benzo[1,2-b:3,4-b']dithiophene (CID
605456) - two structurally distinct isomers (which pair of benzo bonds is
used), both correctly reproduced - plus benzo[1,2-b:4,5-b']difuran (CID
57083408), confirming the mechanism isn't chalcogen-specific.

Scope, deliberately narrow, matching the sibling fusion modules:
- Exactly three rings: one plain benzo (unsubstituted, all-carbon,
  6-membered), plus two *identical* five-membered heteromonocycles (both
  furan, or both thiophene) - a mixed furan+thiophene pair is out of scope
  (needs P-25.3.4.2.1's seniority ordering between different parents, a
  later step).
- Each heteromonocycle is ortho-fused to the benzo ring only (not to each
  other, and not sharing any atom with the other heteromonocycle's own
  fusion bond).
- Each heteromonocycle's own fusion bond must touch its own heteroatom's
  neighbor (the same `_local_numbering`-style restriction the sibling
  modules already carry).
- Completely unsubstituted, no extra atoms, single molecular fragment.

Anything else (different parent rings, 3+ occurrences of a parent, a
parent fused to another parent directly, fusion not reachable from a
heteroatom, any substituent) raises `UnsupportedStructure` and falls
through to other dispatch branches in core.py, exactly like every other
retained/computed-name module in this project.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, ring_cycle
from ._pyridine_heterocycle_fusion import _base_bond_letter, _local_numbering

_PARENT_RING_NAMES = {
    8: "furan",
    16: "thiophene",
}


def _benzo_numbering(graph, benzo_ring, fusion_pairs):
    """Try every (start, direction) numbering of the benzo ring; for each,
    compute both fusion bonds' locant sets and order them (the lower set
    cited first/unprimed, the other primed). Return (numbering,
    first_pair, second_pair) for whichever candidate makes the first-cited
    (unprimed) locant set as low as possible, or None if some candidate's
    fusion atoms aren't a single peripheral bond."""
    cycle = ring_cycle(graph, list(benzo_ring))
    n = len(cycle)
    best = None
    for start in range(n):
        for step in (1, -1):
            numbering = {cycle[(start + step * k) % n]: k + 1 for k in range(n)}
            pair_locants = []
            ok = True
            for pair in fusion_pairs:
                locants = sorted(numbering[atom] for atom in pair)
                if locants[1] - locants[0] != 1 and locants != [1, n]:
                    ok = False
                    break
                pair_locants.append(tuple(locants))
            if not ok:
                continue
            first, second = sorted(pair_locants)
            if best is None or first < best[1]:
                best = (numbering, first, second, pair_locants)
    if best is None:
        return None
    numbering, first, second, pair_locants = best
    first_pair = fusion_pairs[0] if pair_locants[0] == first else fusion_pairs[1]
    second_pair = fusion_pairs[1] if first_pair is fusion_pairs[0] else fusion_pairs[0]
    return numbering, first_pair, second_pair


def find_benzo_bis_heterocycle_fusion_core(mol):
    """Return (benzo_ring, heteroatom_element, [(heteromonocycle_ring,
    heteroatom, fusion_atoms), ...]) if `mol` is exactly one benzo ring
    bridging two identical five-membered heteromonocycles, else None."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 3:
        return None
    sizes = sorted(len(r) for r in atom_rings)
    if sizes != [5, 5, 6]:
        return None
    benzo_ring = next(r for r in atom_rings if len(r) == 6)
    hetero_rings = [r for r in atom_rings if len(r) == 5]

    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic() or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None

    if any(idx for idx in benzo_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6):
        return None

    heteroatoms = []
    for ring in hetero_rings:
        hetero = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
        if len(hetero) != 1 or mol.GetAtomWithIdx(hetero[0]).GetAtomicNum() not in _PARENT_RING_NAMES:
            return None
        heteroatoms.append(hetero[0])
    if mol.GetAtomWithIdx(heteroatoms[0]).GetAtomicNum() != mol.GetAtomWithIdx(heteroatoms[1]).GetAtomicNum():
        return None

    if set(hetero_rings[0]) & set(hetero_rings[1]):
        return None

    fusion_pairs = []
    for ring in hetero_rings:
        shared = set(benzo_ring) & set(ring)
        if len(shared) != 2:
            return None
        a, b = shared
        if b not in {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()}:
            return None
        fusion_pairs.append(shared)
    if fusion_pairs[0] & fusion_pairs[1]:
        return None

    if mol.GetNumAtoms() != len(set(benzo_ring) | set(hetero_rings[0]) | set(hetero_rings[1])):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    element = mol.GetAtomWithIdx(heteroatoms[0]).GetAtomicNum()
    connections = list(zip(hetero_rings, heteroatoms, fusion_pairs))
    return benzo_ring, element, connections


def has_benzo_bis_heterocycle_fusion_name(mol) -> bool:
    return find_benzo_bis_heterocycle_fusion_core(mol) is not None


def name_benzo_bis_heterocycle_fusion(mol) -> str:
    core = find_benzo_bis_heterocycle_fusion_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "this three-ring system is not a supported benzo-bridged "
            "identical-heteromonocycle multiparent fusion (see P-25.3.4.1.3)"
        )
    benzo_ring, element, connections = core
    graph = {atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()}

    fusion_pairs = [pair for _ring, _hetero, pair in connections]
    ordering = _benzo_numbering(graph, benzo_ring, fusion_pairs)
    if ordering is None:
        raise UnsupportedStructure(
            "the two benzo fusion bonds are not both simple peripheral "
            "bonds under a shared numbering (see P-25.3.1.3)"
        )
    benzo_numbering, first_pair, second_pair = ordering
    by_pair = {frozenset(pair): (ring, hetero) for ring, hetero, pair in connections}

    parts = []
    for pair, primed in ((first_pair, False), (second_pair, True)):
        ring, hetero = by_pair[frozenset(pair)]
        parent_numbering = _local_numbering(graph, ring, hetero, pair)
        if parent_numbering is None:
            raise UnsupportedStructure(
                "the fusion bond does not touch a heteromonocycle's own "
                "heteroatom (only the 'b'-lettered fusion is supported, "
                "see P-25.3.1.3)"
            )
        letter_result = _base_bond_letter(graph, ring, hetero, pair)
        if letter_result is None:
            raise UnsupportedStructure(
                "fusion at a heteromonocycle's own heteroatom is not "
                "supported (see P-25.3.1.3)"
            )
        letter, _ = letter_result
        low, high = sorted(pair, key=lambda atom: benzo_numbering[atom])
        prime = "'" if primed else ""
        parts.append(f"{benzo_numbering[low]},{benzo_numbering[high]}-{letter}{prime}")

    parent_name = _PARENT_RING_NAMES[element]
    return f"benzo[{parts[0]}:{parts[1]}]di{parent_name}"
