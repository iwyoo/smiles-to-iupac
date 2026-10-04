"""P-25.3.2.3 orientation and P-25.3.3.1 peripheral numbering of an ortho-fused
tree of five- and six-membered rings, drawn with the permitted shapes.

Every ring is traversed clockwise in one consistent planar embedding. A six-membered ring turns 60 degrees
per edge. A five-membered ring is the house shape (apex up, vertical sides, base edge opposite the apex); its
base is an unfused edge, so every fusion bond of the drawing lies on the hexagonal lattice directions.
Each choice of bases is one drawing; criteria (a)-(d) of P-25.3.2.3.3 pick the orientations, P-25.3.3.1.1 the
start ring and direction, P-25.3.3.1.2 the remaining ties.
"""

import math
from itertools import product

from ._common import UnsupportedStructure
from ._fusion_numbering_general import _HETERO_RANK, _numbering_from, _periphery_cycle_and_fusion_atoms

_EPS = 1e-9
_SHAPE_OFFSETS = {
    5: (0, -90, -150, -210, -270),
    7: (0, -68.7, -90, -150, -210, -270, -291.8),
}
_LATTICE_EDGES = {5: (1, 2, 3, 4), 7: (2, 3, 4, 5)}


def _ring_cycles(mol):
    rings = [list(r) for r in mol.GetRingInfo().AtomRings()]
    if any(len(r) not in (5, 6, 7) for r in rings):
        raise UnsupportedStructure("only five-, six- and seven-membered rings have a measured permitted shape")
    n = len(rings)
    shared = {}
    for i in range(n):
        for j in range(i + 1, n):
            common = set(rings[i]) & set(rings[j])
            if common:
                if len(common) != 2:
                    raise UnsupportedStructure("a peri-fused or spiro ring system has no tree orientation here")
                shared[(i, j)] = shared[(j, i)] = tuple(common)
    if len(shared) // 2 != n - 1:
        raise UnsupportedStructure("a cyclic ring-fusion arrangement has no tree orientation here")
    oriented = {0: rings[0]}
    queue = [0]
    while queue:
        i = queue.pop(0)
        for j in range(n):
            if (i, j) not in shared or j in oriented:
                continue
            u, v = shared[(i, j)]
            cycle = rings[j]
            forward = any(cycle[k] == v and cycle[(k + 1) % len(cycle)] == u for k in range(len(cycle)))
            mine = oriented[i]
            u_to_v = any(mine[k] == u and mine[(k + 1) % len(mine)] == v for k in range(len(mine)))
            wanted_v_to_u_in_j = u_to_v
            oriented[j] = cycle if forward == wanted_v_to_u_in_j else cycle[::-1]
            queue.append(j)
    return [oriented[i] for i in range(n)], shared


def _edge_index(cycle, a, b):
    size = len(cycle)
    for k in range(size):
        if {cycle[k], cycle[(k + 1) % size]} == {a, b}:
            return k
    raise ValueError


def _normals(cycle, parent_edge, parent_normal, base):
    size = len(cycle)
    if size == 6:
        return [(parent_normal - 60 * (e - parent_edge)) % 360 for e in range(6)]
    offsets = _SHAPE_OFFSETS[size]
    shifted = {e: offsets[(e - base) % size] for e in range(size)}
    base_normal = parent_normal - shifted[parent_edge]
    return [(base_normal + shifted[e]) % 360 for e in range(size)]


def _drawings(cycles, shared):
    """Every consistent drawing: ([ring center], [ring normals per edge]) over the choices of five-ring bases."""
    n = len(cycles)
    order = [0]
    parent = {0: None}
    for i in order:
        for j in range(n):
            if (i, j) in shared and j not in parent:
                parent[j] = i
                order.append(j)

    def fused_edges(i):
        return {_edge_index(cycles[i], *shared[(i, j)]) for j in range(n) if (i, j) in shared}

    base_options = []
    for i in range(n):
        size = len(cycles[i])
        if size == 6:
            base_options.append([None])
            continue
        fused = fused_edges(i)
        options = [b for b in range(size) if all((e - b) % size in _LATTICE_EDGES[size] for e in fused)]
        if not options:
            raise UnsupportedStructure("a ring with these fused edges has no permitted-shape placement here")
        base_options.append(options)

    for bases in product(*base_options):
        normals = {}
        center = {0: (0.0, 0.0)}
        for i in order:
            if parent[i] is None:
                edge, normal = 0, 0
                if len(cycles[i]) != 6:
                    edge, normal = bases[i], 270
            else:
                p = parent[i]
                u, v = shared[(i, p)]
                edge = _edge_index(cycles[i], u, v)
                px = _edge_index(cycles[p], u, v)
                normal = (normals[p][px] + 180) % 360
            normals[i] = _normals(cycles[i], edge, normal, bases[i])
            for j in range(n):
                if parent.get(j) == i:
                    u, v = shared[(i, j)]
                    angle = math.radians(normals[i][_edge_index(cycles[i], u, v)])
                    x, y = center[i]
                    center[j] = (x + math.cos(angle), y + math.sin(angle))
        yield [center[i] for i in range(n)], [normals[i] for i in range(n)]


def _row_components(n, shared, cycles, normals, rotation):
    adjacency = {i: set() for i in range(n)}
    for (i, j), (u, v) in shared.items():
        angle = (normals[i][_edge_index(cycles[i], u, v)] - rotation) % 180
        if min(angle, 180 - angle) < _EPS:
            adjacency[i].add(j)
    seen, components = set(), []
    for start in range(n):
        if start in seen:
            continue
        component, stack = [], [start]
        seen.add(start)
        while stack:
            current = stack.pop()
            component.append(current)
            for other in adjacency[current]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        components.append(component)
    return components


def _transform(points, rotation, sx, sy):
    c, s = math.cos(math.radians(-rotation)), math.sin(math.radians(-rotation))
    return [(sx * (x * c - y * s), sy * (x * s + y * c)) for x, y in points]


def _quadrants(points, row):
    ordered = sorted(row, key=lambda i: points[i][0])
    m = len(ordered)
    if m % 2 == 0:
        left, right = points[ordered[m // 2 - 1]], points[ordered[m // 2]]
        cx, cy = (left[0] + right[0]) / 2, left[1]
    else:
        cx, cy = points[ordered[m // 2]]
    upper_right = lower_left = above = 0.0
    for x, y in points:
        x, y = x - cx, y - cy
        on_x, on_y = abs(y) < _EPS, abs(x) < _EPS
        if on_x and on_y:
            upper_right += 0.25
            lower_left += 0.25
            above += 0.5
        elif on_x:
            if x > 0:
                upper_right += 0.5
            else:
                lower_left += 0.5
        elif on_y:
            if y > 0:
                upper_right += 0.5
            else:
                lower_left += 0.5
            above += 0.5 if y > 0 else 0.0
        else:
            if x > 0 and y > 0:
                upper_right += 1.0
            elif x < 0 and y < 0:
                lower_left += 1.0
            if y > 0:
                above += 1.0
    return upper_right, lower_left, above


def _orientations(mol):
    cycles, shared = _ring_cycles(mol)
    n = len(cycles)
    candidates = []
    for centers, normals in _drawings(cycles, shared):
        rotations = sorted({normals[i][_edge_index(cycles[i], *shared[(i, j)])] % 180 for (i, j) in shared})
        for rotation in rotations or [0]:
            components = _row_components(n, shared, cycles, normals, rotation)
            best_row = max(len(c) for c in components)
            for row in (c for c in components if len(c) == best_row):
                for sx, sy in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
                    points = _transform(centers, rotation, sx, sy)
                    ur, ll, above = _quadrants(points, row)
                    candidates.append((best_row, ur, -ll, above, points, sx * sy))
    top = max(c[:4] for c in candidates)
    return cycles, [c for c in candidates if all(abs(a - b) < _EPS for a, b in zip(c[:4], top))]


def oriented_peripheral_numberings(mol, ignore_indicated=False):
    """Every numbering left tied by P-25.3.2.3.3 (a)-(d), P-25.3.3.1.1 and P-25.3.3.1.2, as {atom: locant}."""
    cycles, best = _orientations(mol)
    n = len(cycles)
    cycle, fusion = _periphery_cycle_and_fusion_atoms(mol)
    ring_pairs = {frozenset((c[k], c[(k + 1) % len(c)])): c for c in cycles for k in range(len(c))}
    forward = None
    for k, a in enumerate(cycle):
        b = cycle[(k + 1) % len(cycle)]
        ring = ring_pairs.get(frozenset((a, b)))
        if ring is not None and not (a in fusion and b in fusion):
            idx = ring.index(a)
            forward = ring[(idx + 1) % len(ring)] == b
            break
    clockwise = cycle if forward else cycle[::-1]
    position = {a: i for i, a in enumerate(clockwise)}
    lettered = {a for a in fusion if mol.GetAtomWithIdx(a).GetAtomicNum() == 6}
    numberings = []
    for _row, _ur, _ll, _above, points, handed in best:
        max_y = max(y for _, y in points)
        top_rings = [i for i, (_, y) in enumerate(points) if abs(y - max_y) < _EPS]
        max_x = max(points[i][0] for i in top_rings)
        ring = next(i for i in top_rings if abs(points[i][0] - max_x) < _EPS)
        step = 1 if handed == 1 else -1
        runs = []
        for a in cycles[ring]:
            if a in fusion:
                continue
            before = clockwise[(position[a] - step) % len(clockwise)]
            if before in fusion or before not in cycles[ring]:
                runs.append(position[a])
        if len(runs) != 1:
            raise UnsupportedStructure("the start ring has no single run of nonfusion atoms")
        numberings.append(_numbering_from(clockwise, fusion, runs[0], step, lettered))

    hetero = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 6]
    fusion_carbons = [a for a in fusion if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
    fusion_hetero = [a for a in fusion if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    indicated = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 2]

    def number(text):
        return int(text.rstrip("abcdefgh"))

    def locant_order(text):
        digits = text.rstrip("abcdefgh")
        return int(digits), text[len(digits):]

    def key(locants):
        by_locant = sorted(hetero, key=lambda h: number(locants[h]))
        return (
            [number(locants[h]) for h in by_locant],
            [_HETERO_RANK.get(mol.GetAtomWithIdx(h).GetSymbol(), 99) for h in by_locant],
            sorted(locant_order(locants[a]) for a in fusion_carbons),
            sorted(number(locants[a]) for a in fusion_hetero),
            [] if ignore_indicated else sorted(number(locants[a]) for a in indicated),
        )

    best_key = min(key(x) for x in numberings)
    tied = [x for x in numberings if key(x) == best_key]
    seen, result = set(), []
    for numbering in tied:
        for permutation in mol.GetSubstructMatches(mol, uniquify=False, useChirality=False):
            image = {permutation[a]: loc for a, loc in numbering.items()}
            marker = tuple(sorted(image.items()))
            if marker not in seen:
                seen.add(marker)
                result.append(image)
    return result
