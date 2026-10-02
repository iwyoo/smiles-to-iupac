"""P-25.3.1.3's multi-letter peri-fusion citation (P-25.3.1.1.2's shape,
a third ring spanning three consecutive periphery bonds) for a benzo
ring closing perylene's bay region, its two ends bonding to two new,
mutually-bonded ring atoms. Verified against real PubChem CID 9117
(benzo[ghi]perylene); fully mancude, no indicated-H needed."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._perylene_fusion import _PERYLENE_REF

_CYCLE = [13, 14, 15, 10, 9, 8, 7, 6, 4, 5, 0, 1, 2, 19, 18, 17, 16, 12]
_LETTERS = "abcdefghijklmnopqr"
_N = len(_CYCLE)
_INDEX_OF = {atom: i for i, atom in enumerate(_CYCLE)}


def _find_span(mol):
    matches = mol.GetSubstructMatches(_PERYLENE_REF, useChirality=False, uniquify=False)
    best = None
    for match in matches:
        core = set(match)
        extra = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core]
        if len(extra) != 2:
            continue
        if any(
            mol.GetAtomWithIdx(a).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(a).GetIsAromatic()
            for a in extra
        ):
            continue
        e1, e2 = extra
        if e2 not in [n.GetIdx() for n in mol.GetAtomWithIdx(e1).GetNeighbors()]:
            continue
        target_to_ref = {match[i]: i for i in range(len(match))}
        anchors = []
        for e, other in ((e1, e2), (e2, e1)):
            ring_neighbors = [
                n.GetIdx()
                for n in mol.GetAtomWithIdx(e).GetNeighbors()
                if n.GetIdx() != other and n.GetIdx() in target_to_ref
            ]
            if len(ring_neighbors) != 1:
                anchors = None
                break
            anchors.append(_INDEX_OF[target_to_ref[ring_neighbors[0]]])
        if anchors is None:
            continue
        for start, end in (tuple(anchors), tuple(reversed(anchors))):
            if (end - start) % _N == 3:
                span = "".join(_LETTERS[(start + k) % _N] for k in range(3))
                if best is None or span < best:
                    best = span
    return best


def has_perylene_peri_fusion_name(mol) -> bool:
    return _find_span(mol) is not None


def name_perylene_peri_fusion(mol) -> str:
    span = _find_span(mol)
    if span is None:
        raise UnsupportedStructure(
            "this hexacyclic system is not a supported peri-fused "
            "perylene shape (see P-25.3.1.3)"
        )
    return f"benzo[{span}]perylene"
