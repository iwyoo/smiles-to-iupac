"""Simplified skeletons of phane parent hydrides (P-26.2.1): unbranched acyclic, monocyclic, bicyclic von Baeyer,
polycyclic von Baeyer and monospiro graphs, named 'nonaphane', 'cycloheptaphane', 'bicyclo[6.6.0]tetradecaphane',
'spiro[5.7]tridecaphane' with every numbering the skeletal class allows (P-26.4.1.1)."""

from rdkit import Chem

from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._numerals import numerical_term
from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates
from ._spiro import find_monospiro_atom, iter_monospiro_numberings


def _bare_graph(count, edges):
    skeleton = Chem.RWMol()
    for _ in range(count):
        skeleton.AddAtom(Chem.Atom(0))
    for a, b in edges:
        skeleton.AddBond(a, b, Chem.BondType.SINGLE)
    skeleton = skeleton.GetMol()
    Chem.SanitizeMol(skeleton)
    return skeleton


def _phane_name(parent):
    return parent[:-2] + "phane"


def skeleton_numberings(count, edges):
    """[(name, outer key, node order)] of the skeleton graph with `count` nodes and `edges`; position i of the order
    holds the node with locant i + 1. Empty when no skeletal class of P-26.2.1 describes the graph."""
    neighbours = {i: set() for i in range(count)}
    for a, b in edges:
        neighbours[a].add(b)
        neighbours[b].add(a)
    rings = len(edges) - count + 1
    if rings == 0:
        ends = [i for i, n in neighbours.items() if len(n) == 1]
        if len(ends) != 2 or any(len(n) > 2 for n in neighbours.values()):
            return []
        path = [ends[0]]
        while len(path) < count:
            path.append(next(n for n in neighbours[path[-1]] if n not in path))
        name = numerical_term(count) + "phane"
        return [(name, (), path), (name, (), path[::-1])]
    if any(len(n) < 2 for n in neighbours.values()):
        return []
    if rings == 1:
        if any(len(n) != 2 for n in neighbours.values()):
            return []
        cycle = [0]
        while len(cycle) < count:
            cycle.append(next(n for n in sorted(neighbours[cycle[-1]]) if n not in cycle))
        name = "cyclo" + numerical_term(count) + "phane"
        return [(name, (), (base[start:] + base[:start])) for base in (cycle, cycle[::-1]) for start in range(count)]
    try:
        bare = _bare_graph(count, edges)
    except Exception:
        return []
    if rings == 2:
        spiro = find_monospiro_atom(bare)
        if spiro is not None:
            return [(_phane_name(parent), (), order) for parent, order in iter_monospiro_numberings(bare, spiro)]
        core = find_bicyclic_core(bare)
        if core is None:
            return []
        name = _phane_name(bicyclic_parent_name(core))
        return [(name, (), order) for order in iter_bicyclic_numberings(core)]
    core = find_polycyclic_core(bare, rings)
    if core is None:
        return []
    return [(_phane_name(parent), tuple(outer), order) for order, parent, outer in iter_polycyclic_candidates(core, rings)]
