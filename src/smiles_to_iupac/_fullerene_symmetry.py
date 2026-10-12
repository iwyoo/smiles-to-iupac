"""Point group and ring arrangement of a closed cage of three-coordinate carbon
atoms (P-27.2.1), read from the planar embedding of its 3-connected graph.

A symmetry of the cage is a graph automorphism that keeps the cyclic order of
neighbours (a rotation) or reverses it everywhere (an improper operation); the
Schoenflies symbol follows from the rotation subgroup and the mirror planes."""

from collections import namedtuple

import networkx as nx

Cage = namedtuple("Cage", "symbol surround ambiguous sizes")

_Operation = namedtuple("_Operation", "perm proper")


def _turn(rotation, tail, head, step):
    ring = rotation[tail]
    return ring[(ring.index(head) + step) % len(ring)]


def _edges_of(face):
    return [(face[i], face[(i + 1) % len(face)]) for i in range(len(face))]


def _faces(rotation):
    seen = set()
    faces = []
    for u in rotation:
        for v in rotation[u]:
            face = []
            a, b = u, v
            while (a, b) not in seen:
                seen.add((a, b))
                face.append(a)
                a, b = b, _turn(rotation, b, a, 1)
            if face:
                faces.append(face)
    return faces


def _image(rotation, start, target, proper):
    step = 1 if proper else -1
    dart_image = {start: target}
    vertex_image = {start[0]: target[0], start[1]: target[1]}
    pending = [start]
    while pending:
        a, b = pending.pop()
        a2, b2 = dart_image[(a, b)]
        for dart, image in (
            ((a, _turn(rotation, a, b, 1)), (a2, _turn(rotation, a2, b2, step))),
            ((b, a), (b2, a2)),
        ):
            if dart in dart_image:
                if dart_image[dart] != image:
                    return None
                continue
            dart_image[dart] = image
            pending.append(dart)
            for tail, mapped in zip(dart, image):
                if vertex_image.setdefault(tail, mapped) != mapped:
                    return None
    if len(set(vertex_image.values())) != len(vertex_image):
        return None
    return vertex_image


def _operations(rotation, index):
    vertices = sorted(rotation)
    start = (vertices[0], rotation[vertices[0]][0])
    operations = set()
    for u in vertices:
        for v in rotation[u]:
            for proper in (True, False):
                image = _image(rotation, start, (u, v), proper)
                if image is not None and len(image) == len(vertices):
                    operations.add(_Operation(tuple(index[image[w]] for w in vertices), proper))
    return operations


def _compose(first, second):
    return tuple(second[i] for i in first)


def _order(perm):
    identity = tuple(range(len(perm)))
    power, count = perm, 1
    while power != identity:
        power, count = _compose(power, perm), count + 1
    return count


def _is_mirror(perm, edges, faces):
    return (
        any(perm[i] == i for i in range(len(perm)))
        or any(perm[a] == b and perm[b] == a for a, b in edges)
        or any(frozenset(perm[i] for i in face) == face for face in faces)
    )


def _group_symbol(operations, edges, faces):
    rotations = [op for op in operations if op.proper]
    improper = [op for op in operations if not op.proper]
    size, top = len(rotations), max(_order(op.perm) for op in rotations)
    if size == 1:
        family, n = "C", 1
    elif top == size:
        family, n = "C", size
    elif top * 2 == size:
        family, n = "D", size // 2
    else:
        family, n = {12: "T", 24: "O", 60: "I"}[size], 0
    if not improper:
        return f"{family}{n or ''}"
    mirrors = sum(_order(op.perm) == 2 and _is_mirror(op.perm, edges, faces) for op in improper)
    if family == "C":
        if n == 1:
            return "Cs" if mirrors else "Ci"
        return f"C{n}v" if mirrors == n else f"C{n}h" if mirrors == 1 else f"S{2 * n}"
    if family == "D":
        return f"D{n}h" if mirrors == n + 1 else f"D{n}d"
    if family == "T":
        return "Td" if mirrors == 6 else "Th"
    return f"{family}h"


def _principal_axis(operations, symbol):
    rotations = [op.perm for op in operations if op.proper]
    top = max(_order(p) for p in rotations)
    if top < 3 or symbol[0] in "TOI":
        return None
    generator = next(p for p in rotations if _order(p) == top)
    powers, power = {generator}, generator
    for _ in range(top - 1):
        power = _compose(power, generator)
        powers.add(power)
    if any(_order(p) == top and p not in powers for p in rotations):
        return None
    return generator


def _surround(face_of_dart, face):
    sizes = {len(face_of_dart[(b, a)]) for a, b in _edges_of(face)}
    return sizes.pop() if len(sizes) == 1 else None


def read_cage(graph):
    """Cage(symbol, surround, ambiguous, sizes) of a 3-connected planar cubic `graph`, or None for any other graph.

    `surround` is the size of the rings around the ring on the unique
    principal axis; `ambiguous` is set when that axis meets rings that P-27.2.1
    gives no way to describe by one size."""
    if any(d != 3 for _, d in graph.degree()) or nx.node_connectivity(graph) < 3:
        return None
    planar, embedding = nx.check_planarity(graph)
    if not planar:
        return None
    rotation = {v: list(embedding.neighbors_cw_order(v)) for v in graph}
    index = {v: i for i, v in enumerate(sorted(rotation))}
    faces = _faces(rotation)
    indexed = [frozenset(index[v] for v in f) for f in faces]
    operations = _operations(rotation, index)
    edges = [(index[u], index[v]) for u, v in graph.edges]
    symbol = _group_symbol(operations, edges, indexed)

    surround, ambiguous = None, False
    axis = _principal_axis(operations, symbol)
    if axis is not None:
        face_of_dart = {dart: face for face in faces for dart in _edges_of(face)}
        ends = [f for f, ring in zip(faces, indexed) if frozenset(axis[i] for i in ring) == ring]
        values = {_surround(face_of_dart, f) for f in ends}
        equivalent = len(ends) < 2 or any(
            frozenset(op.perm[i] for i in indexed[faces.index(ends[0])]) == indexed[faces.index(ends[1])]
            for op in operations
        )
        if ends and equivalent and None not in values:
            surround = values.pop()
        elif ends:
            ambiguous = True
    return Cage(symbol, surround, ambiguous, tuple(sorted({len(f) for f in faces})))
