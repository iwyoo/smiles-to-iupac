"""P-25.3.2.3 orientation and P-25.3.3 numbering of an ortho- or ortho- and peri-fused ring system of rings of any
size: rings are laid out from the shapes of P-25.3.2.3.1, rows, quadrants and rings above the row pick the
orientations, numbering starts in the uppermost right-hand ring and runs clockwise, and P-25.3.3.1.2 and
P-25.3.3.2-3 break the remaining ties."""

import math
from collections import defaultdict
from fractions import Fraction
from itertools import product

from rdkit import Chem

from ._common import UnsupportedStructure

HETERO_ORDER = ("F", "Cl", "Br", "I", "O", "S", "Se", "Te", "N", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B", "Al", "Ga", "In", "Tl")
HETERO_RANK = {e: i for i, e in enumerate(HETERO_ORDER)}
_EPS = 1e-6

# Edge normals (turns) of the drawn shapes, clockwise from the right-hand vertical edge (P-25.3.2.3.1)
_HOUSE = {
    5: [Fraction(0), Fraction(3, 4), Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)],
    7: [Fraction(0), Fraction(5, 6), Fraction(2, 3), Fraction(1, 2), Fraction(3, 8), Fraction(1, 4), Fraction(1, 8)],
}


def _shape_tables(n):
    if n in _HOUSE:
        up = _HOUSE[n]
        down = [(-up[(-k) % n]) % 1 for k in range(n)]
        return [up, down]
    return [[Fraction(-k, n) % 1 for k in range(n)]]


class FusedSystem:
    """Rings (clockwise atom cycles in one planar embedding), their shared bonds and the periphery."""

    def __init__(self, mol):
        self.mol = mol
        Chem.GetSymmSSSR(mol)
        rings = [list(r) for r in mol.GetRingInfo().AtomRings()]
        if not rings:
            raise UnsupportedStructure("a ring system without rings has no fusion numbering")
        bond_rings = defaultdict(list)
        for i, ring in enumerate(rings):
            for k in range(len(ring)):
                bond_rings[frozenset((ring[k], ring[(k + 1) % len(ring)]))].append(i)
        if any(len(v) > 2 for v in bond_rings.values()):
            raise UnsupportedStructure("a bond shared by three rings is not a fused ring system")
        self.bond_rings = bond_rings
        self.adj = defaultdict(dict)
        for bond, owners in bond_rings.items():
            if len(owners) == 2:
                i, j = owners
                if j in self.adj[i]:
                    raise UnsupportedStructure("two rings sharing more than one bond are not ortho-fused")
                self.adj[i][j] = bond
                self.adj[j][i] = bond
        self._orient(rings)
        self._periphery()

    def _orient(self, rings):
        n = len(rings)
        done = {0}
        queue = [0]
        while queue:
            i = queue.pop()
            for j, bond in self.adj[i].items():
                a, b = tuple(bond)
                ri, rj = rings[i], rings[j]
                forward_i = ri[(ri.index(a) + 1) % len(ri)] == b
                forward_j = rj[(rj.index(a) + 1) % len(rj)] == b
                if j in done:
                    if forward_i == forward_j:
                        raise UnsupportedStructure("an inconsistent ring embedding (not a planar fused system)")
                    continue
                if forward_i == forward_j:
                    rings[j] = rj[::-1]
                done.add(j)
                queue.append(j)
        if len(done) != n:
            raise UnsupportedStructure("a disconnected ring system")
        self.rings = rings
        self.ring_sets = [set(r) for r in rings]
        self.tree = []
        seen = {0}
        queue = [0]
        while queue:
            i = queue.pop(0)
            for j in sorted(self.adj[i]):
                if j not in seen:
                    seen.add(j)
                    self.tree.append((i, j))
                    queue.append(j)

    def _periphery(self):
        successor = {}
        for i, ring in enumerate(self.rings):
            for k in range(len(ring)):
                a, b = ring[k], ring[(k + 1) % len(ring)]
                if len(self.bond_rings[frozenset((a, b))]) == 1:
                    if a in successor:
                        raise UnsupportedStructure("a ring system whose periphery touches itself")
                    successor[a] = b
        if not successor:
            raise UnsupportedStructure("a closed ring surface has no periphery")
        start = min(successor)
        cycle = [start]
        while successor[cycle[-1]] != start:
            cycle.append(successor[cycle[-1]])
        if len(cycle) != len(successor):
            raise UnsupportedStructure("a ring system with more than one periphery")
        self.periphery = cycle
        counts = defaultdict(int)
        for ring in self.rings:
            for a in ring:
                counts[a] += 1
        self.fusion = {a for a in cycle if counts[a] > 1}
        self.interior = sorted(a for a in counts if a not in set(cycle))

    def edge_index(self, ring, bond):
        cyc = self.rings[ring]
        n = len(cyc)
        return next(k for k in range(n) if frozenset((cyc[k], cyc[(k + 1) % n])) == bond)


def _layouts(system):
    """Every drawing of the ring system: per ring the shape role choices, giving ring centres and the normal of
    each ring edge. Yields (centres, normals)."""
    rings = system.rings
    n = len(rings)
    parent = {j: i for i, j in system.tree}
    order = [0] + [j for _, j in system.tree]
    choices = []
    for r in order:
        size = len(rings[r])
        if size in _HOUSE:
            tables = _shape_tables(size)
            choices.append([(t, j) for t in range(len(tables)) for j in range(size)])
        else:
            choices.append([(0, 0)])
    seen = set()
    found = []
    for pick in product(*choices):
        centres = {0: (0.0, 0.0)}
        normals = {}
        for idx, r in enumerate(order):
            size = len(rings[r])
            variant, role = pick[idx]
            table = _shape_tables(size)[variant]
            if r == 0:
                p_edge = 0
                p_normal = Fraction(0)
            else:
                p = parent[r]
                bond = system.adj[r][p]
                p_edge = system.edge_index(r, bond)
                p_normal = (normals[(p, system.edge_index(p, bond))] + Fraction(1, 2)) % 1
            for e in range(size):
                normals[(r, e)] = (p_normal + table[(e - p_edge + role) % size] - table[role]) % 1
            for child in system.adj[r]:
                if parent.get(child) == r:
                    bond = system.adj[r][child]
                    e = system.edge_index(r, bond)
                    theta = 2 * math.pi * float(normals[(r, e)])
                    dist = 1.0
                    cx, cy = centres[r]
                    centres[child] = (cx + dist * math.cos(theta), cy + dist * math.sin(theta))
        key = tuple((round(centres[r][0], 4), round(centres[r][1], 4)) for r in range(n)) + tuple(
            sorted((k, v) for k, v in normals.items())
        )
        if key in seen:
            continue
        seen.add(key)
        consistent = all(
            (normals[(i, system.edge_index(i, bond))] - normals[(j, system.edge_index(j, bond))] - Fraction(1, 2)) % 1 == 0
            for i in system.adj
            for j, bond in system.adj[i].items()
        )
        found.append((consistent, centres, normals))
    drawn = [d for d in found if d[0]] or found
    for _, centres, normals in drawn:
        yield centres, normals


def _weights(delta):
    if abs(delta) < _EPS:
        return 0.5, 0.5
    return (1.0, 0.0) if delta > 0 else (0.0, 1.0)


def _row_components(system, normals, residue):
    link = defaultdict(set)
    for i in system.adj:
        for j, bond in system.adj[i].items():
            if (normals[(i, system.edge_index(i, bond))] - residue) % Fraction(1, 2) == 0:
                link[i].add(j)
    seen = set()
    rows = []
    for start in range(len(system.rings)):
        if start in seen:
            continue
        comp, stack = [], [start]
        seen.add(start)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in link[cur]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        rows.append(comp)
    return rows


def _orientations(system, with_key=False):
    """Winning drawings of P-25.3.2.3.3 (a)-(d): list of (ring centres, orientation preserving); the winning key
    (negated row length, negated upper-right quadrant rings, lower-left rings, negated rings above) on request."""
    results = []
    best = None
    for centres, normals in _layouts(system):
        residues = {
            (normals[(i, system.edge_index(i, bond))]) % Fraction(1, 2) for i in system.adj for bond in system.adj[i].values()
        } or {Fraction(0)}
        for residue in residues:
            theta = -2 * math.pi * float(residue)
            rotated = {
                r: (c[0] * math.cos(theta) - c[1] * math.sin(theta), c[0] * math.sin(theta) + c[1] * math.cos(theta))
                for r, c in centres.items()
            }
            for row in _row_components(system, normals, residue):
                for sx, sy in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
                    placed = {r: (sx * x, sy * y) for r, (x, y) in rotated.items()}
                    members = sorted(row, key=lambda r: placed[r][0])
                    m = len(members)
                    if m % 2:
                        cx, cy = placed[members[m // 2]]
                    else:
                        left, right = placed[members[m // 2 - 1]], placed[members[m // 2]]
                        cx, cy = (left[0] + right[0]) / 2, left[1]
                    upper_right = lower_left = above = 0.0
                    for x, y in placed.values():
                        up, down = _weights(y - cy)
                        rt, lt = _weights(x - cx)
                        upper_right += rt * up
                        lower_left += lt * down
                        above += up
                    key = (-m, -round(upper_right, 6), round(lower_left, 6), -round(above, 6))
                    if best is None or key < best:
                        best, results = key, []
                    if key == best:
                        results.append((placed, sx * sy > 0))
    return (results, best) if with_key else results


def orientation_key(mol):
    """P-25.3.2.3.3 key of the preferred orientation; the smaller key is the senior orientation."""
    return _orientations(FusedSystem(mol), with_key=True)[1]


def _locant_key(text):
    digits = ""
    for ch in text:
        if ch.isdigit():
            digits += ch
        else:
            break
    rest = text[len(digits):]
    letters = "".join(c for c in rest if c.isalpha())
    sup = "".join(c for c in rest if c.isdigit())
    return int(digits), letters, int(sup) if sup else 0


def _number_periphery(system, walk, start_index):
    mol = system.mol
    locants = {}
    counter = 0
    letter = 0
    previous_fusion = False
    for step in range(len(walk)):
        atom = walk[(start_index + step) % len(walk)]
        carbon = mol.GetAtomWithIdx(atom).GetAtomicNum() == 6
        if atom in system.fusion and carbon:
            letter = letter + 1 if previous_fusion else 1
            locants[atom] = f"{counter}{chr(96 + letter)}"
            previous_fusion = True
        else:
            counter += 1
            locants[atom] = str(counter)
            previous_fusion = False
    return locants, counter


def _number_interior(system, locants, counter):
    mol = system.mol
    carbons = {a for a in system.interior if mol.GetAtomWithIdx(a).GetAtomicNum() == 6}
    hetero = [a for a in system.interior if a not in carbons]
    adjacency = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()}
    peripheral = set(locants)

    def nearest(atom):
        frontier, seen, steps = {atom}, {atom}, 0
        while frontier:
            steps += 1
            reached = {nb for cur in frontier for nb in adjacency[cur] if nb in peripheral}
            if reached:
                return steps, min(reached, key=lambda x: _locant_key(locants[x]))
            frontier = {nb for cur in frontier for nb in adjacency[cur] if nb not in seen}
            seen |= frontier
        raise UnsupportedStructure("an interior atom not connected to the periphery")

    for atom in sorted(hetero, key=lambda a: (nearest(a)[0], _locant_key(locants[nearest(a)[1]]), HETERO_RANK.get(mol.GetAtomWithIdx(a).GetSymbol(), 99))):
        counter += 1
        locants[atom] = str(counter)
    for atom in carbons:
        steps, origin = nearest(atom)
        locants[atom] = f"{locants[origin]}{steps}"


def _distance(system, source, targets):
    adjacency = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in system.mol.GetAtoms()}
    frontier, seen, steps = {source}, {source}, 0
    while frontier:
        if frontier & targets:
            return steps
        steps += 1
        frontier = {nb for cur in frontier for nb in adjacency[cur] if nb not in seen}
        seen |= frontier
    return steps


def _start_atoms(system, ordered, ring):
    atoms = system.ring_sets[ring]
    n = len(ordered)
    free = lambda a: a not in system.fusion
    starts = [a for i, a in enumerate(ordered) if a in atoms and free(a) and not (ordered[i - 1] in atoms and free(ordered[i - 1]))]
    if starts:
        return starts
    arc = next(i for i, a in enumerate(ordered) if a in atoms and ordered[i - 1] not in atoms)
    for step in range(n):
        atom = ordered[(arc + step) % n]
        if free(atom):
            return [atom]
    raise UnsupportedStructure("a fused ring system without a nonfusion atom")


def _helicene_ends(system):
    """The two terminal rings of an unbranched chain of six or more hexagons that all turn the same way, else None."""
    n = len(system.rings)
    if n < 6 or any(len(r) != 6 for r in system.rings):
        return None
    if any(len(system.adj[i]) > 2 for i in range(n)) or sum(len(system.adj[i]) for i in range(n)) != 2 * (n - 1):
        return None
    ends = [i for i in range(n) if len(system.adj[i]) == 1]
    order, previous = [ends[0]], None
    while len(order) < n:
        following = [j for j in system.adj[order[-1]] if j != previous][0]
        previous = order[-1]
        order.append(following)
    turns = set()
    for before, ring, after in zip(order, order[1:], order[2:]):
        a = system.edge_index(ring, system.adj[ring][before])
        b = system.edge_index(ring, system.adj[ring][after])
        turns.add((b - a) % 6)
    return ends if len(turns) == 1 and turns <= {2, 4} else None


def fused_numberings(mol):
    """Every numbering {atom: locant text} of the fused ring system left tied after P-25.3.2.3.3, P-25.3.3.1 and
    P-25.3.3.1.2 (a)-(e) and P-25.3.3.3.2; indicated hydrogen (f) is left to the caller."""
    system = FusedSystem(mol)
    produced = {}
    ends = _helicene_ends(system)
    if ends is not None:
        # P-25.3.3.1.1: a terminal ring is drawn upper right and numbering begins in it
        drawings = [(None, preserving) for preserving in (True, False)]
    else:
        drawings = _orientations(system)
    for placed, preserving in drawings:
        ordered = system.periphery if preserving else system.periphery[::-1]
        if ends is not None:
            index_of = {a: i for i, a in enumerate(ordered)}
            for ring in ends:
                for atom in _start_atoms(system, ordered, ring):
                    locants, counter = _number_periphery(system, ordered, index_of[atom])
                    _number_interior(system, locants, counter)
                    produced[tuple(sorted(locants.items()))] = locants
            continue
        index_of = {a: i for i, a in enumerate(ordered)}
        top = max(y for _, y in placed.values())
        uppermost = [r for r, (x, y) in placed.items() if abs(y - top) < _EPS]
        right = max(placed[r][0] for r in uppermost)
        for ring in [r for r in uppermost if abs(placed[r][0] - right) < _EPS]:
            for atom in _start_atoms(system, ordered, ring):
                locants, counter = _number_periphery(system, ordered, index_of[atom])
                _number_interior(system, locants, counter)
                produced[tuple(sorted(locants.items()))] = locants
    if not produced:
        raise UnsupportedStructure("no numbering of this fused ring system")
    hetero = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 6]
    fusion_carbons = [a for a in system.fusion if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
    fusion_hetero = [a for a in hetero if a in system.fusion]
    interior_hetero = [a for a in hetero if a in system.interior]
    interior_carbon = [a for a in system.interior if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]

    def key(locants):
        by_locant = sorted(hetero, key=lambda h: _locant_key(locants[h]))
        lowest_fusion = min(fusion_carbons, key=lambda a: _locant_key(locants[a]), default=None)
        return (
            [_locant_key(locants[h]) for h in by_locant],
            [HETERO_RANK.get(mol.GetAtomWithIdx(h).GetSymbol(), 99) for h in by_locant],
            sorted(_locant_key(locants[a]) for a in fusion_carbons),
            sorted(_locant_key(locants[a]) for a in fusion_hetero),
            sorted(_distance(system, h, {lowest_fusion}) for h in interior_hetero) if lowest_fusion is not None else [],
            sorted(_locant_key(locants[a]) for a in interior_carbon),
        )

    ranked = list(produced.values())
    best = min(key(c) for c in ranked)
    tied = [c for c in ranked if key(c) == best]
    final, seen = [], set()
    for numbering in tied:
        for perm in mol.GetSubstructMatches(mol, uniquify=False, useChirality=False):
            image = {perm[a]: loc for a, loc in numbering.items()}
            ident = tuple(sorted(image.items()))
            if ident not in seen:
                seen.add(ident)
                final.append(image)
    return final
