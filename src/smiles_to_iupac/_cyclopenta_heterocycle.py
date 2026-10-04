"""Cyclopentane-ring-fused five-membered heteromonocycles named cyclopenta[b|c] + component
(P-25.3.1.3): the fusion letter is b when a fusion atom is bonded to the heteroatom, else c.
Whole-system numbering and indicated hydrogen follow `_benzo_heterocycle`.
"""

from ._benzo_heterocycle import _HETEROATOMS, _component, _mancude_requirement, _saturated
from ._common import UnsupportedStructure
from ._fusion_numbering_general import _HETERO_RANK, general_peripheral_numberings


def _match(mol):
    if mol.GetRingInfo().NumRings() != 2 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if len(rings[0]) != 5 or len(rings[1]) != 5 or len(rings[0] & rings[1]) != 2:
        return None
    if mol.GetNumAtoms() != len(rings[0] | rings[1]):
        return None
    carbocycles = [r for r in rings if all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in r)]
    if len(carbocycles) != 1:
        return None
    other = rings[0] if carbocycles[0] is rings[1] else rings[1]
    hetero = [a for a in other if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(hetero) != 1 or mol.GetAtomWithIdx(hetero[0]).GetSymbol() not in _HETEROATOMS:
        return None
    sat = _saturated(mol, rings[0] | rings[1])
    if len(sat) != _mancude_requirement(mol, rings[0] | rings[1]):
        return None
    return carbocycles[0], other, hetero[0], sat


def has_cyclopenta_heterocycle_name(mol) -> bool:
    return _match(mol) is not None


def name_cyclopenta_heterocycle(mol) -> str:
    carbocycle, other, hetero, sat = _match(mol)
    numberings = general_peripheral_numberings(mol, ignore_indicated=True)
    if not numberings:
        raise UnsupportedStructure("this cyclopenta-fused ring has no supported numbering")
    locant = lambda numbering, a: int(numbering[a].rstrip("abcdefgh"))
    numbering = min(numberings, key=lambda n: (_HETERO_RANK[mol.GetAtomWithIdx(hetero).GetSymbol()], locant(n, hetero), sorted(locant(n, a) for a in sat)))
    ih = sorted(locant(numbering, a) for a in sat)
    fusion = other & carbocycle
    letter = "b" if any(n.GetIdx() == hetero for a in fusion for n in mol.GetAtomWithIdx(a).GetNeighbors()) else "c"
    ih_text = ",".join(f"{p}H" for p in ih) + "-" if ih else ""
    return f"{ih_text}cyclopenta[{letter}]{_component(mol, other)}"
