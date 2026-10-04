"""P-25.3.3.1 orientation and peripheral numbering of ortho-fused systems made only of six-membered rings,
drawn on a hexagonal lattice: rings in the horizontal row, rings in the upper right and lower left quadrants and
rings above the row decide the orientation, numbering starts in the upper right ring and runs clockwise, then
heteroatoms, fusion carbon atoms and indicated hydrogen take the lowest locants."""

import math
from collections import defaultdict

from ._fusion_numbering_general import _HETERO_RANK

_ROOT3 = math.sqrt(3)
_EPS = 1e-6


def _ring_cycle(mol, atoms):
    adjacency = defaultdict(list)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in atoms and b in atoms:
            adjacency[a].append(b)
            adjacency[b].append(a)
    start = min(atoms)
    cycle, previous, current = [start], None, start
    while True:
        following = next(x for x in adjacency[current] if x != previous)
        if following == start:
            return cycle
        cycle.append(following)
        previous, current = current, following


def _vertex(center, index):
    angle = math.radians(90 - 60 * index)
    return center[0] + math.cos(angle), center[1] + math.sin(angle)


def _embed(mol):
    """({atom: (x, y)}, [ring centres], [ring atom sets]) of the fused system on a hexagonal lattice, or None."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if len(rings) < 2 or any(len(r) != 6 for r in rings):
        return None
    if any(sum(a in r for r in rings) > 3 for a in range(mol.GetNumAtoms())):
        return None
    cycles = [_ring_cycle(mol, r) for r in rings]
    centers = {0: (0.0, 0.0)}
    clockwise = {0: cycles[0]}
    coordinates = {atom: _vertex((0.0, 0.0), i) for i, atom in enumerate(cycles[0])}
    queue = [0]
    while queue:
        current = queue.pop()
        for other, ring in enumerate(rings):
            if other == current or other in centers or len(rings[current] & ring) != 2:
                continue
            shared = rings[current] & ring
            order = clockwise[current]
            k = next(i for i in range(6) if order[i] in shared and order[(i + 1) % 6] in shared)
            angle = math.radians(60 - 60 * k)
            center = (centers[current][0] + _ROOT3 * math.cos(angle), centers[current][1] + _ROOT3 * math.sin(angle))
            centers[other] = center
            a, b = order[k], order[(k + 1) % 6]

            def vertex_index(atom):
                dx, dy = coordinates[atom][0] - center[0], coordinates[atom][1] - center[1]
                return round((90 - math.degrees(math.atan2(dy, dx))) / 60) % 6

            cycle = cycles[other]
            ia, ib = cycle.index(a), cycle.index(b)
            step = 1 if (ib - ia) % 6 == 1 else -1
            if vertex_index(b) != (vertex_index(a) + 1) % 6:
                step = -step
            ordered = [None] * 6
            for t in range(6):
                ordered[(vertex_index(a) + t) % 6] = cycle[(ia + step * t) % 6]
            clockwise[other] = ordered
            for i, atom in enumerate(ordered):
                coordinates.setdefault(atom, _vertex(center, i))
            queue.append(other)
    if len(centers) != len(rings):
        return None
    return coordinates, [centers[i] for i in range(len(rings))], rings


def _transform(point, rotation, mirror):
    x, y = point
    if mirror:
        x = -x
    angle = math.radians(60 * rotation)
    return x * math.cos(angle) - y * math.sin(angle), x * math.sin(angle) + y * math.cos(angle)


def _weights(delta):
    if abs(delta) < _EPS:
        return 0.5, 0.5
    return (1.0, 0.0) if delta > 0 else (0.0, 1.0)


def _rows(centers):
    by_height = defaultdict(list)
    for i, (x, y) in enumerate(centers):
        by_height[round(y / 1.5)].append((x, i))
    rows = []
    for members in by_height.values():
        members.sort()
        run = [members[0]]
        for item in members[1:]:
            if abs(item[0] - run[-1][0] - _ROOT3) < _EPS:
                run.append(item)
            else:
                rows.append([i for _, i in run])
                run = [item]
        rows.append([i for _, i in run])
    return rows


def _periphery(mol, atoms=None):
    bond_rings = defaultdict(int)
    for ring in mol.GetRingInfo().BondRings():
        for bond in ring:
            bond_rings[bond] += 1
    adjacency = defaultdict(list)
    for bond in mol.GetBonds():
        if bond_rings[bond.GetIdx()] == 1:
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            adjacency[a].append(b)
            adjacency[b].append(a)
    atoms = sorted(adjacency)
    fusion = {a for a in atoms if sum(1 for r in mol.GetRingInfo().AtomRings() if a in r) > 1}
    start = min(atoms)
    cycle, previous, current = [start], None, start
    while True:
        following = next(x for x in adjacency[current] if x != previous)
        if following == start:
            break
        cycle.append(following)
        previous, current = current, following
    return cycle, fusion


def _is_helicene(centers):
    """A chain of angularly fused rings that all turn the same way (the helicene orientation rule differs)."""
    neighbours = {
        i: [j for j in range(len(centers)) if j != i and abs(math.dist(centers[i], centers[j]) - _ROOT3) < _EPS]
        for i in range(len(centers))
    }
    if any(len(v) > 2 for v in neighbours.values()):
        return False
    ends = [i for i, v in neighbours.items() if len(v) == 1]
    if len(ends) != 2:
        return False
    order, previous = [ends[0]], None
    while len(order) < len(centers):
        following = [j for j in neighbours[order[-1]] if j != previous][0]
        previous = order[-1]
        order.append(following)
    turns = []
    for a, b, c in zip(order, order[1:], order[2:]):
        v1 = (centers[b][0] - centers[a][0], centers[b][1] - centers[a][1])
        v2 = (centers[c][0] - centers[b][0], centers[c][1] - centers[b][1])
        cross = v1[0] * v2[1] - v1[1] * v2[0]
        if abs(cross) < _EPS:
            return False
        turns.append(cross > 0)
    return len(set(turns)) == 1


def _number_interior(mol, locants, interior):
    """P-25.3.3.3: an interior atom takes the locant of the nearest peripheral atom plus the bond count as a
    superscript (3a1), the lowest peripheral locant when several are equally near."""
    if not interior:
        return
    pending = set(interior)
    frontier = set(locants)
    steps = 0
    while pending:
        steps += 1
        reached = {}
        for atom in frontier:
            for n in mol.GetAtomWithIdx(atom).GetNeighbors():
                idx = n.GetIdx()
                if idx in pending:
                    reached.setdefault(idx, []).append(atom)
        if not reached:
            return
        for idx, sources in reached.items():
            origin = min(sources, key=lambda a: _locant_key(locants[a]))
            locants[idx] = f"{locants[origin]}{steps}"
            pending.discard(idx)
        frontier = set(reached)


def _locant_key(text):
    digits = ""
    for ch in text:
        if ch.isdigit():
            digits += ch
        else:
            break
    return int(digits), text[len(digits):]


def _letters(count):
    return chr(ord("a") + count - 1)


def hex_numberings(mol, ignore_indicated=True):
    """Every tied numbering {atom: locant text} of an ortho-fused all-six-membered-ring system, or None."""
    embedding = _embed(mol)
    if embedding is None:
        return None
    coordinates, centers, rings = embedding
    if len(rings) >= 5 and _is_helicene(centers):
        return None
    candidates = []
    best_key = None
    for mirror in (False, True):
        for rotation in range(6):
            placed = [_transform(c, rotation, mirror) for c in centers]
            for row in _rows(placed):
                cx = (placed[row[0]][0] + placed[row[-1]][0]) / 2
                cy = placed[row[0]][1]
                upper_right = lower_left = above = 0.0
                for x, y in placed:
                    up, down = _weights(y - cy)
                    right, left = _weights(x - cx)
                    upper_right += right * up
                    lower_left += left * down
                    above += up
                key = (-len(row), -upper_right, lower_left, -above)
                if best_key is None or key < best_key:
                    best_key, candidates = key, []
                if key == best_key:
                    candidates.append((mirror, rotation))
    cycle, fusion = _periphery(mol)
    interior = [a for a in range(mol.GetNumAtoms()) if a not in set(cycle)]
    hetero = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 6]
    fusion_carbons = [a for a in fusion if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
    indicated = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 2]
    results = {}
    for mirror, rotation in candidates:
        placed = {a: _transform(p, rotation, mirror) for a, p in coordinates.items()}
        area = sum(
            placed[cycle[i]][0] * placed[cycle[(i + 1) % len(cycle)]][1]
            - placed[cycle[(i + 1) % len(cycle)]][0] * placed[cycle[i]][1]
            for i in range(len(cycle))
        )
        walk = cycle if area < 0 else cycle[::-1]
        ring_height = {}
        for ring_index, ring in enumerate(rings):
            ring_height[ring_index] = _transform(centers[ring_index], rotation, mirror)
        order = sorted(range(len(rings)), key=lambda r: (-round(ring_height[r][1] / 1.5), -ring_height[r][0]))
        start_atom = None
        for ring_index in order:
            free = [a for a in rings[ring_index] if a not in fusion]
            if not free:
                continue
            n = len(walk)
            for i, atom in enumerate(walk):
                if atom in free and walk[(i - 1) % n] not in free:
                    start_atom = atom
                    break
            if start_atom is not None:
                break
        if start_atom is None:
            continue
        start = walk.index(start_atom)
        locants, counter, letter = {}, 0, 0
        previous_fusion = False
        for step in range(len(walk)):
            atom = walk[(start + step) % len(walk)]
            if atom in fusion:
                letter = letter + 1 if previous_fusion else 1
                locants[atom] = f"{counter}{_letters(letter)}"
                previous_fusion = True
            else:
                counter += 1
                locants[atom] = str(counter)
                previous_fusion = False
        _number_interior(mol, locants, interior)
        results[tuple(sorted(locants.items()))] = locants

    def locant_order(text):
        return _locant_key(text)

    def sort_key(locants):
        hetero_key = sorted(locant_order(locants[h]) for h in hetero)
        rank_key = [
            _HETERO_RANK.get(mol.GetAtomWithIdx(h).GetSymbol(), 99)
            for h in sorted(hetero, key=lambda h: locant_order(locants[h]))
        ]
        fusion_key = sorted(locant_order(locants[a]) for a in fusion_carbons)
        indicated_key = [] if ignore_indicated else sorted(locant_order(locants[h]) for h in indicated)
        return hetero_key, rank_key, fusion_key, indicated_key

    if not results:
        return None
    ranked = list(results.values())
    best = min(sort_key(c) for c in ranked)
    tied = [c for c in ranked if sort_key(c) == best]
    seen, final = set(), []
    for numbering in tied:
        for permutation in mol.GetSubstructMatches(mol, uniquify=False, useChirality=False):
            image = {permutation[a]: loc for a, loc in numbering.items()}
            key = tuple(sorted(image.items()))
            if key not in seen:
                seen.add(key)
                final.append(image)
    return final
