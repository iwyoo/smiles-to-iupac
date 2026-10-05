"""Systematic atom numbering of fullerene cages (P-27.3; Fu-3.1 of the IUPAC fullerene recommendations).

Fu-3.1 never defines a "contiguous spiral pathway" (it defers to CAS, ref. 3 there, and says its rules suffice only
for C60-Ih and C70-D5h(6)). Here it is the atom walk that, on arriving at an atom from its predecessor, takes the
first unvisited neighbour in a fixed rotational sense, and is contiguous when it visits every cage atom.
"""

from itertools import combinations

import networkx as nx
from rdkit import Chem

_CLASSES = ("ring", "bond", "atom")


class _Cage:
    def __init__(self, graph):
        _, embedding = nx.check_planarity(graph)
        self.graph = graph
        self.rotation = {v: list(embedding.neighbors_cw_order(v)) for v in graph}
        self.faces = _faces(embedding)
        self.face_sizes_at = {v: sorted((len(f) for f in self.faces if v in f), reverse=True) for v in graph}
        self.rotations = _rotations(self.rotation)
        self.features = _features(self)


def _faces(embedding):
    seen = set()
    faces = []
    for u in embedding:
        for v in embedding[u]:
            if (u, v) not in seen:
                faces.append(embedding.traverse_face(u, v, mark_half_edges=seen))
    return faces


def _rotations(rotation):
    """Orientation-preserving automorphisms of the cage as atom permutations (identity included)."""
    first = next(iter(rotation))
    anchor = rotation[first][0]
    found = []
    for x in rotation:
        for y in rotation[x]:
            image = _propagate(rotation, first, anchor, x, y)
            if image is not None:
                found.append(image)
    return found


def _propagate(rotation, first, anchor, x, y):
    image = {first: x}
    stack = [(first, anchor, x, y)]
    while stack:
        a, b, x, y = stack.pop()
        row_a, row_x = rotation[a], rotation[x]
        if len(row_a) != len(row_x):
            return None
        shift = row_x.index(y) - row_a.index(b)
        for i, n in enumerate(row_a):
            target = row_x[(i + shift) % len(row_x)]
            if n in image:
                if image[n] != target:
                    return None
            else:
                image[n] = target
                stack.append((n, a, target, x))
    if len(set(image.values())) != len(rotation):
        return None
    return image


def _stabilizer(cage, atoms):
    return [g for g in cage.rotations if frozenset(g[a] for a in atoms) == frozenset(atoms)]


def _features(cage):
    """(class, order, atoms) of every ring, bond and atom lying on a proper rotation axis."""
    out = []
    for face in cage.faces:
        order = len(_stabilizer(cage, face))
        if order > 1:
            out.append(("ring", order, tuple(face)))
    for u, v in cage.graph.edges:
        order = len(_stabilizer(cage, (u, v)))
        if order > 1:
            out.append(("bond", order, (u, v)))
    for atom in cage.graph:
        order = len(_stabilizer(cage, (atom,)))
        if order > 1:
            out.append(("atom", order, (atom,)))
    return out


def _spiral(cage, start, second, clockwise):
    path = [start, second]
    visited = set(path)
    step = 1 if clockwise else -1
    while len(path) < cage.graph.number_of_nodes():
        row = cage.rotation[path[-1]]
        index = row.index(path[-2])
        following = [row[(index + step * k) % len(row)] for k in range(1, len(row))]
        unvisited = [n for n in following if n not in visited]
        if not unvisited:
            return None
        path.append(unvisited[0])
        visited.add(unvisited[0])
    return path


def _starts(cage, feature):
    kind, _, atoms = feature
    if kind == "bond":
        u, v = atoms
        return [(u, v), (v, u)]
    if kind == "atom":
        (a,) = atoms
        return [(a, n) for n in cage.rotation[a]]
    ring = set(atoms)
    return [(a, n) for a in atoms for n in cage.graph[a] if n in ring]


def _winds_ring(cage, ring, start, second, clockwise):
    row = cage.rotation[second]
    step = 1 if clockwise else -1
    return row[(row.index(start) + step) % len(row)] in ring


def _rank(cage, feature):
    kind, _, atoms = feature
    signatures = sorted((tuple(cage.face_sizes_at[a]) for a in atoms), reverse=True)
    return (len(atoms) if kind == "ring" else 0, tuple(signatures))


def _antipode(cage, feature):
    """Atoms of the feature at the other end of the rotation axis through `feature`."""
    turn = next(g for g in _stabilizer(cage, feature[2]) if any(g[a] != a for a in cage.graph))
    atoms = set(feature[2])
    for other in cage.features:
        if set(other[2]) != atoms and frozenset(turn[a] for a in other[2]) == frozenset(other[2]):
            return set(other[2])
    return set()


def _paths(cage, feature):
    kind = feature[0]
    ring = set(feature[2])
    paths = []
    for start, second in _starts(cage, feature):
        for clockwise in (True, False):
            if kind == "ring" and not _winds_ring(cage, ring, start, second, clockwise):
                continue
            path = _spiral(cage, start, second, clockwise)
            if path is not None:
                paths.append(path)
    return paths


def cage_numberings(mol, pool, attach):
    """`Numbering`s of a C60-Ih/C70-D5h(6) cage as a parent or a group with free valences (P-29.3.4.1), its
    saturated atoms as added hydrogen and hydro prefixes, e.g. '(C60-Ih)[5,6]fulleren-1(9H)-yl'."""
    from ._common import UnsupportedStructure
    from ._fullerene import _cage_graph, numbered_cage_stem
    from ._ring_diyl_numbering import Numbering

    stem = numbered_cage_stem(mol, pool)
    if stem is None:  # Fu-3.1 gives rules only for C60-Ih and C70-D5h(6); the CAS rules for other cages are external
        raise UnsupportedStructure(
            "fullerene locants cannot be derived: Fu-3.1 states its rules suffice only for the C60-Ih and "
            "C70-D5h(6) cages and that more rules are needed for other fullerenes"
        )
    free = set(attach)
    saturated = {
        a for a in pool
        if mol.GetAtomWithIdx(a).GetTotalNumHs() > 0
        or mol.GetAtomWithIdx(a).GetFormalCharge() != 0
        or mol.GetAtomWithIdx(a).GetNumRadicalElectrons() > 0
        or any(n.GetIdx() not in pool for n in mol.GetAtomWithIdx(a).GetNeighbors())
    } | free
    others = sorted(saturated - free)
    # P-58.2: one added hydrogen when the free-valence atoms are odd in number, hydro prefixes for the other pairs
    added_count = len(free) % 2
    if (len(others) - added_count) % 2 or len(others) < added_count:
        raise UnsupportedStructure("the saturated atoms of this fullerene cannot be paired as hydro/added hydrogen")

    numberings = fullerene_numberings(_cage_graph(mol, pool))
    if not numberings:
        raise UnsupportedStructure("no contiguous spiral pathway starts on a rotation axis of this fullerene")
    out = []
    for position_of in numberings:
        for added_atoms in combinations(others, added_count):
            added = tuple(sorted(position_of[a] for a in added_atoms))
            hydro = tuple(sorted(position_of[a] for a in others if a not in added_atoms))
            out.append(Numbering(position_of, _cage_text(stem, added, hydro), unsat_key=(added, hydro)))
    return out


_PREFIX_ONLY_ELEMENTS = {6, 9, 17, 35, 53}


def name_cage_parent(mol):
    """Name of a C60-Ih/C70-D5h(6) cage carrying only hydrocarbon and halogen substituents, as the parent
    hydride with hydro prefixes, e.g. '1-(trifluoromethyl)-1,9-dihydro(C60-Ih)[5,6]fullerene' (P-6, P-31.1.4)."""
    from ._common import UnsupportedStructure, adjacency
    from ._diester_ring_diyl import _system_of, evaluate_skeleton
    from ._fullerene import is_fullerene_cage

    pentagon = next(r for r in mol.GetRingInfo().AtomRings() if len(r) == 5)
    rings, cage = _system_of(mol, pentagon[0])
    if not is_fullerene_cage(mol, cage):
        raise UnsupportedStructure("this ring system is not a fullerene cage")
    if any(
        mol.GetAtomWithIdx(a).GetAtomicNum() not in _PREFIX_ONLY_ELEMENTS
        or mol.GetAtomWithIdx(a).GetFormalCharge()
        or mol.GetAtomWithIdx(a).GetNumRadicalElectrons()
        for a in range(mol.GetNumAtoms())
        if a not in cage
    ) or len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure(
            "a fullerene parent is named here only with hydrocarbon and halogen substituents; characteristic "
            "groups on the cage need suffix naming (P-6)"
        )
    found = evaluate_skeleton(mol, adjacency(mol), "ring", rings, cage, [], set(), "")
    if found is None:
        raise UnsupportedStructure("this substituted fullerene has no supported name yet")
    return found[1]


def _cage_text(stem, added, hydro):
    from ._ring_diyl_numbering import _hydro_text, _tail_added

    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        rest = _tail_added(stem, locants, valence, suffix, added)
        prefix = _hydro_text(hydro)
        return prefix + ("-" if prefix and rest[0].isdigit() else "") + rest

    return text


def fullerene_numberings(graph):
    """Every systematic numbering (atom -> locant) of the cage `graph` that ties under Fu-3.1; None when no
    contiguous spiral pathway starts on a rotation axis."""
    cage = _Cage(graph)
    for order in sorted({f[1] for f in cage.features}, reverse=True):
        for kind in _CLASSES:
            found = [
                (f, p) for f in cage.features if f[1] == order and f[0] == kind for p in [_paths(cage, f)] if p
            ]
            if not found:
                continue
            best = max(_rank(cage, f) for f, _ in found)
            scored = []
            for feature, paths in found:
                if _rank(cage, feature) == best:
                    antipode = _antipode(cage, feature)
                    for path in paths:
                        distance = min((nx.shortest_path_length(cage.graph, path[-1], a) for a in antipode), default=0)
                        scored.append((distance, tuple(path)))
            closest = min(d for d, _ in scored)
            return [{atom: i + 1 for i, atom in enumerate(p)} for p in sorted({p for d, p in scored if d == closest})]
    return None
