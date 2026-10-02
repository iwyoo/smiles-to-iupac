"""P-25.3.1.3's multi-letter peri-fusion citation applied twice on the
same base: two benzo rings closing both bay regions of chrysene, at
P-25.3.1.3's 'def' and 'mno' spans (chrysene's own C2 symmetry maps one
span onto the other). Verified by InChI match against real PubChem CID
9118 (anthanthrene, PIN dibenzo[def,mno]chrysene)."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._chrysene_fusion import _CHRYSENE_REF

_CYCLE = [5, 0, 1, 2, 3, 9, 10, 11, 12, 13, 14, 15, 16, 17, 8, 7, 6, 4]
_LETTERS = "abcdefghijklmnopqr"
_N = len(_CYCLE)
_INDEX_OF = {atom: i for i, atom in enumerate(_CYCLE)}


def _bridge_span(mol, core, target_to_ref, e1, e2):
    if e2 not in [n.GetIdx() for n in mol.GetAtomWithIdx(e1).GetNeighbors()]:
        return None
    anchors = []
    for e, other in ((e1, e2), (e2, e1)):
        ring_neighbors = [
            n.GetIdx()
            for n in mol.GetAtomWithIdx(e).GetNeighbors()
            if n.GetIdx() != other and n.GetIdx() in target_to_ref
        ]
        if len(ring_neighbors) != 1:
            return None
        anchors.append(_INDEX_OF[target_to_ref[ring_neighbors[0]]])
    for start, end in (tuple(anchors), tuple(reversed(anchors))):
        if (end - start) % _N == 3:
            return "".join(_LETTERS[(start + k) % _N] for k in range(3))
    return None


def _find_spans(mol):
    matches = mol.GetSubstructMatches(_CHRYSENE_REF, useChirality=False, uniquify=False)
    best = None
    for match in matches:
        core = set(match)
        extra = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core]
        if len(extra) != 4:
            continue
        if any(
            mol.GetAtomWithIdx(a).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(a).GetIsAromatic()
            for a in extra
        ):
            continue
        target_to_ref = {match[i]: i for i in range(len(match))}
        unpaired = set(extra)
        spans = []
        ok = True
        while unpaired:
            e1 = next(iter(unpaired))
            partners = [
                n.GetIdx()
                for n in mol.GetAtomWithIdx(e1).GetNeighbors()
                if n.GetIdx() in unpaired and n.GetIdx() != e1
            ]
            if len(partners) != 1:
                ok = False
                break
            e2 = partners[0]
            span = _bridge_span(mol, core, target_to_ref, e1, e2)
            if span is None:
                ok = False
                break
            spans.append(span)
            unpaired -= {e1, e2}
        if not ok or len(spans) != 2:
            continue
        citation = ",".join(sorted(spans))
        if best is None or citation < best:
            best = citation
    return best


def has_anthanthrene_fusion_name(mol) -> bool:
    return _find_spans(mol) is not None


def name_anthanthrene_fusion(mol) -> str:
    spans = _find_spans(mol)
    if spans is None:
        raise UnsupportedStructure(
            "this system is not a supported double peri-fused chrysene "
            "shape (see P-25.3.1.3)"
        )
    return f"dibenzo[{spans}]chrysene"
