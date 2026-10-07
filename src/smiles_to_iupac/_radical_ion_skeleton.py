"""Radical ions on one saturated unbranched hydrocarbon skeleton (P-75.2): the anionic or cationic suffixes follow the
parent hydride, then the radical suffixes; radical centres take the lowest locants."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._numerals import alkane_name, multiplying_prefix


def _centres(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    ions = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if (
        not radicals
        or not ions
        or any(a.GetNumRadicalElectrons() != 1 or a.GetFormalCharge() for a in radicals)
        or len({a.GetFormalCharge() for a in ions}) != 1
        or abs(ions[0].GetFormalCharge()) != 1
        or any(a.GetNumRadicalElectrons() for a in ions)
    ):
        return None
    return radicals, ions


def _shape(mol):
    if (
        len(Chem.GetMolFrags(mol)) != 1
        or mol.GetNumAtoms() < 2
        or any(a.GetAtomicNum() != 6 or a.GetIsotope() or a.GetIsAromatic() or a.GetDegree() > 2 for a in mol.GetAtoms())
        or any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds())
    ):
        return None
    rings = mol.GetRingInfo().NumRings()
    if rings > 1 or (rings == 1 and mol.GetBonds().__len__() != mol.GetNumAtoms()):
        return None
    return rings


def has_skeleton_radical_ion_shape(mol) -> bool:
    return _centres(mol) is not None and _shape(mol) is not None


def _orders(mol, ring):
    graph = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()}
    size = mol.GetNumAtoms()
    if not ring:
        start = next(i for i, n in graph.items() if len(n) == 1)
        order, previous = [start], None
        while len(order) < size:
            nxt = [n for n in graph[order[-1]] if n != previous]
            previous = order[-1]
            order.append(nxt[0])
        return [order, order[::-1]]
    cycle = [next(iter(graph))]
    previous = None
    while len(cycle) < size:
        nxt = [n for n in graph[cycle[-1]] if n != previous]
        previous = cycle[-1]
        cycle.append(nxt[0])
    return [[cycle[(s + step * k) % size] for k in range(size)] for s in range(size) for step in (1, -1)]


def _suffix(count, singular):
    return singular if count == 1 else multiplying_prefix(count) + singular


def name_skeleton_radical_ion(mol) -> str:
    centres = _centres(mol)
    ring = _shape(mol)
    if centres is None or ring is None:
        raise UnsupportedStructure("this radical ion is not a plain unbranched hydrocarbon skeleton")
    radicals, ions = centres
    if any(a.GetIdx() == r.GetIdx() for a in ions for r in radicals):
        raise UnsupportedStructure("a radical and an ionic centre on one atom are not supported yet")
    best = None
    for order in _orders(mol, bool(ring)):
        locant = {atom: i + 1 for i, atom in enumerate(order)}
        key = (sorted(locant[a.GetIdx()] for a in radicals), sorted(locant[a.GetIdx()] for a in ions))
        if best is None or key < best:
            best = key
    radical_locants, ion_locants = best
    ending = "id" if ions[0].GetFormalCharge() < 0 else "ium"
    stem = ("cyclo" if ring else "") + alkane_name(mol.GetNumAtoms())
    ion_text, radical_text = _suffix(len(ions), ending), _suffix(len(radicals), "yl")
    if ion_text[0] not in "aeiou":
        stem += ""
    else:
        stem = stem[:-1]
    return f"{stem}-{','.join(map(str, ion_locants))}-{ion_text}-{','.join(map(str, radical_locants))}-{radical_text}"
