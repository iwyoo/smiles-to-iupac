"""General fullerene-isomer identification via the Fowler-Manolopoulos ring
spiral algorithm (P-27.3's "Numbering of fullerenes" references this as the
basis of the systematic point-group-qualified names used throughout P-27.2.3,
e.g. the already-shipped `(C70-D5h(6))[5,6]fullerene`).

`_fullerene.py`'s exact-SMILES table only covers cage sizes small enough to
have a single (or a single *isolated*) IPR isomer (C60, C70, C76). Beyond
that, C84 alone has 24 distinct isolated-pentagon-rule isomers, and knowing
which one a given structure is requires comparing its own topology against a
real reference atlas -- exactly the "isomer-atlas/spiral-code cross-reference
tool" that issue #997 said didn't exist in this project. This module is that
tool, built from the published algorithm rather than any single lookup
shortcut:

1. `_planar_faces`: get the fullerene's 12 pentagons + N hexagons as an exact
   combinatorial fact, not a heuristic. A fullerene cage graph is 3-connected
   and planar, so by Whitney's uniqueness theorem it has exactly one
   embedding up to reflection -- `networkx.check_planarity` finds it, and
   `PlanarEmbedding.traverse_face` walks each face's boundary directly from
   that embedding's rotation system. This is deliberately *not* done via
   RDKit's ring perception: RDKit's symmetrized SSSR recovers the full face
   set for highly-symmetric cages (as used in `_fullerene.py`'s C60/C70/C76
   checks) but can under-count faces for the low point-group-symmetry
   isomers (C1, C2, Cs, ...) that dominate C84's 24-isomer list, since
   symmetrization only adds the *extra* rings a real symmetry axis demands.
2. `_canonical_spiral`: a spiral is the face sequence you'd get by "peeling"
   the polyhedron in a single continuous strip from a chosen starting face
   and edge. Extracting one from an *existing* structure (as opposed to
   generating a fullerene from a candidate spiral, the more commonly
   documented direction) is a breadth-first search of the planar dual graph
   (faces as nodes, edge-sharing faces adjacent) where each face's neighbors
   are visited starting immediately after the parent edge and continuing in
   a fixed rotational direction -- i.e. the BFS order is fully determined by
   (starting face, starting edge, handedness). Every (face, edge, handedness)
   combination gives *a* valid spiral; the canonical one is whichever gives
   the lexicographically smallest sequence of face sizes (pentagons sort
   before hexagons), which is the standard atlas canonicalization. This
   project's own implementation of that BFS reproduces, byte for byte, the
   published canonical spirals for all three cages `_fullerene.py` already
   ships by independent means:
   - C60-Ih: pentagons at faces (1,7,9,11,13,15,18,20,22,24,26,32) of 32
   - C70-D5h(6): (1,7,9,11,13,15,27,29,31,33,35,37) of 37
   - C76-D2: (1,7,9,11,13,18,26,31,33,35,37,39) of 40
   (cross-checked against the CTCP "Program Fullerene" IPR database,
   https://ctcp.massey.ac.nz/ -- the same Fowler-Manolopoulos-derived dataset
   cited below for C84.) A real, independently-fetched PubChem structure
   (CID 133108900, C84, InChIKey FQRWAZOLUJHNDT-UHFFFAOYSA-N) spirals to
   (1,7,10,13,19,22,25,28,30,34,37,44), which is exactly isomer #24 (D6h) in
   that same database -- confirmed again by a real paper's own title, "A
   minor isomer of C84 fullerene, D6h-C84(24)" (ScienceDirect), naming that
   same isomer the same way.
3. `_C84_IPR_ISOMERS`: the 24 IPR isomers' canonical spirals and point
   groups, transcribed from the CTCP database (`c084IPR.database`), cross-
   checked row for row against the independent MSU nanotube fullerene
   database (nanotube.msu.edu/fullerene). The literature convention (seen in
   the D6h-C84(24) paper title above, and in CAS-style names generally) is
   to always cite the parenthesized global canonical index, even for a
   point-group symbol that happens to be unique among the 24 -- unlike this
   project's existing bare `(C76-D2)[5,6]fullerene` (C76 has only one IPR
   isomer total in the relevant sense, so the question never arises there).

Out of scope, same as `_fullerene.py`'s own exact-match table: non-IPR
isomers (chemically unstable, essentially never the structure a real input
SMILES encodes), any cage size other than 60/70/76/84 (no reference atlas
data embedded here for them -- see this module's docstring in the project's
Issue tracker for how to extend it), and any substituted/modified fullerene
at any size (P-27.3's full atom-by-atom numbering system, needed only for
locants on substituents or modifications, is separate unimplemented work).
"""

from collections import deque

import networkx as nx


def _planar_faces(mol):
    """Return the list of faces (each a list of atom indices) of `mol`'s
    bond graph, or None if the graph isn't planar (not a fullerene cage)."""
    graph = nx.Graph()
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    is_planar, embedding = nx.check_planarity(graph)
    if not is_planar:
        return None
    faces = []
    visited_half_edges = set()
    for u in embedding:
        for v in embedding[u]:
            if (u, v) in visited_half_edges:
                continue
            faces.append(embedding.traverse_face(u, v, mark_half_edges=visited_half_edges))
    return faces


def _face_adjacency(faces):
    """For each face (by index), its neighboring faces in the same cyclic
    order as its own boundary edges."""
    half_edge_owner = {}
    for i, face in enumerate(faces):
        n = len(face)
        for k in range(n):
            u, v = face[k], face[(k + 1) % n]
            half_edge_owner[(u, v)] = i
    neighbors = []
    for face in faces:
        n = len(face)
        row = []
        for k in range(n):
            u, v = face[k], face[(k + 1) % n]
            row.append(half_edge_owner[(v, u)])
        neighbors.append(row)
    return neighbors


def _spiral_from(neighbors, start_face, start_rotation, mirror):
    """BFS the dual graph from `start_face`, starting at boundary-edge offset
    `start_rotation` and using the reversed (mirror) rotational direction if
    `mirror`. Returns the face-visit order, or None if some face is
    unreachable (not expected for a connected fullerene cage)."""
    total = len(neighbors)
    rows = {i: (list(reversed(row)) if mirror else row) for i, row in enumerate(neighbors)}
    order = [start_face]
    visited = {start_face}
    parent_of = {}
    queue = deque()

    start_row = rows[start_face]
    rotated = start_row[start_rotation:] + start_row[:start_rotation]
    for neighbor in rotated:
        if neighbor not in visited:
            visited.add(neighbor)
            order.append(neighbor)
            parent_of[neighbor] = start_face
            queue.append(neighbor)

    while queue:
        face = queue.popleft()
        row = rows[face]
        parent_index = row.index(parent_of[face])
        rotated = row[parent_index + 1:] + row[:parent_index + 1]
        for neighbor in rotated:
            if neighbor not in visited:
                visited.add(neighbor)
                order.append(neighbor)
                parent_of[neighbor] = face
                queue.append(neighbor)

    return order if len(order) == total else None


def _canonical_spiral(faces):
    """The lexicographically smallest face-size sequence reachable by a
    spiral BFS from any (starting face, starting edge, handedness)."""
    neighbors = _face_adjacency(faces)
    sizes = [len(face) for face in faces]
    best = None
    for start_face in range(len(faces)):
        for start_rotation in range(len(neighbors[start_face])):
            for mirror in (False, True):
                order = _spiral_from(neighbors, start_face, start_rotation, mirror)
                if order is None:
                    continue
                code = tuple(sizes[i] for i in order)
                if best is None or code < best:
                    best = code
    return best


def _pentagon_positions(spiral_code):
    return tuple(i + 1 for i, size in enumerate(spiral_code) if size == 5)


# The 24 isolated-pentagon-rule isomers of C84, transcribed from the CTCP
# "Program Fullerene" IPR database (c084IPR.database), in the database's own
# canonical numbering (identical to Fowler & Manolopoulos's atlas numbering
# per that program's own documentation). Keyed by canonical pentagon
# positions (1-indexed, among the cage's 44 faces) for direct lookup against
# `_pentagon_positions(_canonical_spiral(...))`.
_C84_IPR_ISOMERS = {
    (1, 7, 9, 11, 13, 18, 24, 35, 38, 40, 42, 44): (1, "D2"),
    (1, 7, 9, 11, 13, 24, 28, 30, 36, 40, 42, 44): (2, "C2"),
    (1, 7, 9, 11, 14, 22, 27, 29, 31, 35, 41, 43): (3, "Cs"),
    (1, 7, 9, 11, 14, 22, 27, 30, 35, 39, 41, 43): (4, "D2d"),
    (1, 7, 9, 11, 14, 23, 28, 30, 36, 40, 42, 44): (5, "D2"),
    (1, 7, 9, 11, 22, 24, 26, 28, 30, 32, 36, 44): (6, "C2v"),
    (1, 7, 9, 11, 22, 24, 26, 28, 30, 32, 42, 44): (7, "C2v"),
    (1, 7, 9, 12, 14, 20, 26, 29, 33, 37, 40, 42): (8, "C2"),
    (1, 7, 9, 12, 14, 20, 27, 29, 32, 35, 41, 43): (9, "C2"),
    (1, 7, 9, 12, 14, 20, 27, 29, 33, 35, 40, 43): (10, "Cs"),
    (1, 7, 9, 12, 20, 24, 26, 28, 30, 33, 36, 44): (11, "C2"),
    (1, 7, 9, 12, 20, 24, 26, 28, 30, 33, 42, 44): (12, "C1"),
    (1, 7, 9, 12, 20, 24, 26, 28, 30, 34, 41, 44): (13, "C2"),
    (1, 7, 9, 12, 20, 24, 26, 28, 33, 36, 39, 41): (14, "Cs"),
    (1, 7, 9, 12, 21, 24, 26, 28, 30, 32, 42, 44): (15, "Cs"),
    (1, 7, 9, 13, 20, 22, 25, 28, 30, 34, 37, 44): (16, "Cs"),
    (1, 7, 9, 13, 20, 22, 25, 28, 30, 37, 41, 43): (17, "C2v"),
    (1, 7, 9, 13, 20, 22, 25, 28, 34, 37, 39, 41): (18, "C2v"),
    (1, 7, 9, 13, 20, 22, 26, 28, 30, 34, 36, 44): (19, "D3d"),
    (1, 7, 9, 13, 20, 23, 25, 28, 33, 37, 39, 41): (20, "Td"),
    (1, 7, 10, 12, 14, 18, 26, 31, 33, 37, 39, 42): (21, "D2"),
    (1, 7, 10, 13, 18, 22, 25, 27, 31, 34, 38, 44): (22, "D2"),
    (1, 7, 10, 13, 18, 22, 25, 27, 31, 38, 41, 43): (23, "D2d"),
    (1, 7, 10, 13, 19, 22, 25, 28, 30, 34, 37, 44): (24, "D6h"),
}


def match_c84_isomer(mol):
    """If `mol` is an all-carbon, degree-3, 84-atom cage whose canonical
    spiral matches one of C84's 24 known IPR isomers, return the
    `"(C84-<point group>(<index>))[5,6]fullerene"` PIN fragment. Otherwise
    None (not a C84 cage, not planar/3-connected, or an isomer -- IPR or
    not -- outside this table)."""
    if mol.GetNumAtoms() != 84:
        return None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C" or atom.GetDegree() != 3:
            return None

    faces = _planar_faces(mol)
    if faces is None or len(faces) != 44:
        return None
    sizes = [len(face) for face in faces]
    if sizes.count(5) != 12 or sizes.count(6) != 32:
        return None

    spiral = _canonical_spiral(faces)
    pentagons = _pentagon_positions(spiral)
    entry = _C84_IPR_ISOMERS.get(pentagons)
    if entry is None:
        return None
    index, point_group = entry
    return f"(C84-{point_group}({index}))[5,6]fullerene"
