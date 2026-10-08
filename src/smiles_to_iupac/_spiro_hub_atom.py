"""Three monocyclic rings through one spiro atom of nonstandard bonding number (P-24.8.1.3): 'spiro' before the name of
the alicyclic system with the same number of atoms, the spiro atom's locant as a superscript each time it is revisited,
and the heteroatoms, the spiro atom among them, as 'a' prefixes: 1,4,6,9,10,13-hexaoxa-5λ6-thiaspiro[4.4^5.4^5]tridecane."""

from itertools import permutations, product

from rdkit import Chem

from ._common import UnsupportedStructure, ring_cycle
from ._numerals import alkane_name, numerical_term

_REPLACEMENT = {8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza"}
_SENIORITY = (8, 16, 34, 52, 7)
_HUB = {16: 6, 34: 6, 52: 6}


def _shape(mol):
    """(hub atom, [ring cycles]) for a molecule made of three saturated monocycles sharing one atom, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetRingInfo().NumRings() != 3:
        return None
    rings = [list(r) for r in mol.GetRingInfo().AtomRings()]
    common = set(rings[0]) & set(rings[1]) & set(rings[2])
    if len(common) != 1 or any(len(set(a) & set(b)) != 1 for a, b in ((rings[0], rings[1]), (rings[0], rings[2]), (rings[1], rings[2]))):
        return None
    hub = mol.GetAtomWithIdx(next(iter(common)))
    if hub.GetAtomicNum() not in _HUB or hub.GetDegree() != 6 or hub.GetFormalCharge() or hub.GetTotalNumHs():
        return None
    if sum(len(r) for r in rings) - 2 != mol.GetNumAtoms():
        return None
    for atom in mol.GetAtoms():
        if atom.GetIdx() == hub.GetIdx():
            continue
        if atom.GetAtomicNum() not in {6, *_REPLACEMENT} or atom.GetDegree() != 2 or atom.GetFormalCharge() or atom.GetIsotope():
            return None
    if any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()):
        return None
    adjacency = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()}
    cycles = []
    for ring in rings:
        inside = {a: [n for n in adjacency[a] if n in ring] for a in ring}
        cycles.append(ring_cycle(inside, ring))
    return hub.GetIdx(), cycles


def has_spiro_hub_atom_shape(mol) -> bool:
    return _shape(mol) is not None


def _arms(cycle, hub):
    start = cycle.index(hub)
    ordered = cycle[start + 1 :] + cycle[:start]
    return ordered, ordered[::-1]


def name_spiro_hub_atom(mol) -> str:
    found = _shape(mol)
    if found is None:
        raise UnsupportedStructure("not three monocycles through one nonstandard spiro atom")
    hub, cycles = found
    best = None
    for order in permutations(cycles):
        sizes = [len(c) - 1 for c in order]
        if sizes != sorted(sizes):
            continue
        for directions in product((0, 1), repeat=3):
            locant, number = {}, 1
            for cycle, direction in zip(order, directions):
                for atom in _arms(cycle, hub)[direction]:
                    locant[atom] = number
                    number += 1
                if cycle is order[0]:
                    locant[hub] = number
                    number += 1
            hetero = {a: mol.GetAtomWithIdx(a).GetAtomicNum() for a in locant if mol.GetAtomWithIdx(a).GetAtomicNum() != 6}
            key = (
                sorted(locant[a] for a in hetero),
                [tuple(sorted(locant[a] for a, z in hetero.items() if z == element)) for element in _SENIORITY],
            )
            if best is None or key < best[0]:
                best = (key, locant, sizes, hetero)
    _, locant, sizes, hetero = best
    pieces = []
    for element in _SENIORITY:
        atoms = sorted((a for a, z in hetero.items() if z == element), key=locant.get)
        if not atoms:
            continue
        cited = ",".join(f"{locant[a]}λ{_HUB[element]}" if a == hub else str(locant[a]) for a in atoms)
        multiplier = numerical_term(len(atoms)) if len(atoms) > 1 else ""
        pieces.append(f"{cited}-{multiplier}{_REPLACEMENT[element]}")
    spiro = locant[hub]
    descriptor = f"[{sizes[0]}.{sizes[1]}^{spiro}.{sizes[2]}^{spiro}]"
    return f"{'-'.join(pieces)}spiro{descriptor}{alkane_name(mol.GetNumAtoms())}"
