"""Polyspiro systems of monocycles joined through spiro atoms in any tree arrangement (P-24.2.2, P-24.2.3, P-24.2.4.1):
dispiro, trispiro ... nonaspiro with the von Baeyer spiro descriptor, skeletal replacement prefixes in front.

The descriptor is read on a closed walk that starts in a terminal ring and returns to the first spiro atom: a ring is
entered at a spiro atom, run round with each branching spiro atom followed by an excursion into the ring attached to it,
and left at the entry atom again. Numbers count the atoms of each stretch; a stretch ending on an atom met before cites
that atom's locant as a superscript. Lowest locants go to the spiro atoms, then the descriptor numbers in order of
citation, then the heteroatoms."""

from itertools import product

from rdkit import Chem

from ._common import adjacency, ring_cycle
from ._numerals import alkane_name, numerical_term
from ._vb_ring_assembly import _PREFIX, _SENIORITY, _VALENCE


_ORDER = sorted(_SENIORITY, key=_SENIORITY.get)


def _tree(mol):
    """(rings as cycles, {spiro atom: (ring i, ring j)}) of a molecule that is only monocycles joined by spiro atoms."""
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumAtoms() != sum(1 for a in mol.GetAtoms() if a.IsInRing()):
        return None
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if len(rings) < 3 or any(b.GetBondTypeAsDouble() != 1.0 or b.GetIsAromatic() for b in mol.GetBonds()):
        return None
    shared = {}
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            common = rings[i] & rings[j]
            if len(common) > 1:
                return None
            if common:
                (atom,) = common
                shared.setdefault(atom, set()).update((i, j))
    if any(len(pair) != 2 for pair in shared.values()) or len(shared) != len(rings) - 1:
        return None
    degree = {i: sum(i in pair for pair in shared.values()) for i in range(len(rings))}
    seen, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for pair in shared.values():
            if i in pair:
                for j in pair - {i}:
                    if j not in seen:
                        seen.add(j)
                        stack.append(j)
    if len(seen) != len(rings) or any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z != 6 and (z not in _PREFIX or atom.GetTotalValence() != _VALENCE[z]):
            return None
    graph = adjacency(mol)
    cycles = [ring_cycle(graph, list(ring)) for ring in rings]
    return cycles, {atom: tuple(sorted(pair)) for atom, pair in shared.items()}, degree


def has_spiro_tree_shape(mol) -> bool:
    return _tree(mol) is not None


def _ring_walks(cycle, entry, stops):
    """Both directions round `cycle` from `entry`: (atoms after the entry in order, back at the entry)."""
    start = cycle.index(entry)
    forward = cycle[start + 1:] + cycle[:start]
    return [forward, forward[::-1]]


def _walk(cycles, spiro, ring, entry, parent):
    """Yield event lists for a ring run round from `entry`; event = ('plain', atom) | ('stop', atom, sub-events)."""
    options = []
    for atoms in _ring_walks(cycles[ring], entry, spiro):
        pieces = []
        for atom in atoms:
            child = None
            if atom in spiro:
                (other,) = [r for r in spiro[atom] if r != ring]
                if other != parent:
                    child = other
            if child is None:
                pieces.append([[("plain", atom)]])
            else:
                pieces.append([[("stop", atom, sub)] for sub in _walk(cycles, spiro, child, atom, ring)])
        for combination in product(*pieces):
            options.append([event for events in combination for event in events])
    return options


def _flatten(events, entry, numbering, descriptor, spiro_locants):
    """Number the atoms met on a walk and emit descriptor items; a stretch ends at a stop or back at `entry`."""
    count = 0
    for event in events:
        if event[0] == "plain":
            count += 1
            numbering.append(event[1])
        else:
            atom, sub = event[1], event[2]
            numbering.append(atom)
            spiro_locants.append(len(numbering))
            descriptor.append((count, None))
            count = 0
            _flatten(sub, atom, numbering, descriptor, spiro_locants)
    descriptor.append((count, entry))


def _candidates(cycles, spiro, degree):
    terminals = [i for i, d in degree.items() if d == 1]
    for start_ring in terminals:
        (entry,) = [a for a, pair in spiro.items() if start_ring in pair]
        (second,) = [r for r in spiro[entry] if r != start_ring]
        start_cycle = cycles[start_ring]
        index = start_cycle.index(entry)
        around = start_cycle[index + 1:] + start_cycle[:index]
        for first_atoms in (around, around[::-1]):
            for sub in _walk(cycles, spiro, second, entry, start_ring):
                numbering = list(first_atoms)
                descriptor = [(len(first_atoms), None)]
                numbering.append(entry)
                spiro_locants = [len(numbering)]
                tail = []
                _flatten(sub, entry, numbering, tail, spiro_locants)
                descriptor.extend(tail)
                yield numbering, descriptor, spiro_locants


def name_spiro_tree(mol) -> str:
    cycles, spiro, degree = _tree(mol)
    best = None
    for numbering, descriptor, spiro_locants in _candidates(cycles, spiro, degree):
        position = {atom: i + 1 for i, atom in enumerate(numbering)}
        items = []
        for count, back in descriptor:
            items.append(f"{count}^{position[back]}" if back is not None else str(count))
        hetero = [(mol.GetAtomWithIdx(a).GetAtomicNum(), position[a]) for a in numbering if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
        locant_set = tuple(sorted(p for _, p in hetero))
        by_element = tuple(tuple(sorted(p for z, p in hetero if z == element)) for element in _ORDER)
        key = (tuple(sorted(spiro_locants)), tuple(c for c, _ in descriptor), tuple(spiro_locants), locant_set, by_element)
        if best is None or key < best[0]:
            best = (key, items, hetero)
    _, items, hetero = best
    hydrocarbon = f"{numerical_term(len(spiro))}spiro[{'.'.join(items)}]{alkane_name(mol.GetNumAtoms())}"
    if not hetero:
        return hydrocarbon
    pieces = []
    for z in sorted({z for z, _ in hetero}, key=_SENIORITY.get):
        locants = sorted(p for element, p in hetero if element == z)
        multiplier = numerical_term(len(locants)) if len(locants) > 1 else ""
        pieces.append(f"{','.join(map(str, locants))}-{multiplier}{_PREFIX[z]}")
    return "-".join(pieces) + hydrocarbon
