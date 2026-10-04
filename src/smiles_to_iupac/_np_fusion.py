"""Mancude rings fused to a natural-product parent: fusion prefixes with primed locants (P-101.5.1)."""

from dataclasses import dataclass

from ._common import UnsupportedStructure
from ._np_core import FACES, loc_key

_PRIORITY = ["O", "S", "Se", "Te", "N", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B"]
_HETERO_PREFIX = {
    "O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha", "As": "arsa", "Sb": "stiba",
    "Bi": "bisma", "Si": "sila", "Ge": "germa", "Sn": "stanna", "Pb": "plumba", "B": "bora",
}
_CARBO = {3: "cyclopropa", 4: "cyclobuta", 5: "cyclopenta", 6: "benzo", 7: "cyclohepta", 8: "cycloocta"}
_STEM = {3: "irene", 4: "ete", 5: "ole", 6: "ine", 7: "epine", 8: "ocine"}
_RETAINED = {
    (5, ("O",)): "furo", (5, ("S",)): "thieno", (5, ("Se",)): "selenopheno", (5, ("Te",)): "telluropheno",
    (5, ("N",)): "pyrrolo", (6, ("N",)): "pyrido", (6, ("O",)): "pyrano", (6, ("S",)): "thiopyrano",
    (5, ("N", "N")): "pyrazolo", (6, ("N", "N")): "pyrazino",
}
_PRIMES = ["′", "″", "‴"]


@dataclass
class Fused:
    prefix: str
    name_key: str
    atoms: set
    shared: list
    indicated: list


def _hw_name(size, hetero):
    """Hantzsch-Widman fusion prefix of a mancude ring and whether its hetero locants are cited."""
    elements = tuple(sorted(hetero, key=_PRIORITY.index))
    if (size, elements) in _RETAINED:
        return _RETAINED[(size, elements)], False
    if not elements:
        return _CARBO[size], False
    counts = {}
    for el in elements:
        counts[el] = counts.get(el, 0) + 1
    prefix = ""
    mult = {2: "di", 3: "tri", 4: "tetra"}
    for el in sorted(counts, key=_PRIORITY.index):
        prefix += (mult.get(counts[el], "") if counts[el] > 1 else "") + _HETERO_PREFIX[el]
    stem = _STEM.get(size)
    if stem is None:
        raise UnsupportedStructure("a ring size without a Hantzsch-Widman stem")
    if size == 3 and elements == ("N",):
        text = "azirine"
    elif size == 6 and elements[0] in ("N",) and all(e == "N" for e in elements):
        raise UnsupportedStructure("a polyaza six-membered component without a retained name")
    else:
        glue = prefix
        if glue.endswith("a") and stem[0] in "aeiou":
            glue = glue[:-1]
        text = glue + stem
        text = text.replace("aa", "a")
    name = text[:-1] + "o" if text.endswith("e") else text + "o"
    return name, len(elements) > 1 or (len(elements) == 1 and False)


def _cycle_from(atoms, adj, start):
    cycle, previous = [start], None
    while True:
        nxt = [n for n in adj[cycle[-1]] if n in atoms and n != previous and n not in cycle[1:]]
        if not nxt or (len(cycle) > 1 and nxt[0] == start):
            break
        previous = cycle[-1]
        cycle.append(nxt[0])
        if len(cycle) > len(atoms):
            break
    return cycle


def ring_cycle(atoms, adj):
    atoms = set(atoms)
    start = next(iter(sorted(atoms)))
    cycle = [start]
    previous = None
    while True:
        options = [n for n in adj[cycle[-1]] if n in atoms and n != previous]
        if len(cycle) == len(atoms):
            break
        step = next((n for n in options if n not in cycle), None)
        if step is None:
            return None
        previous = cycle[-1]
        cycle.append(step)
    return cycle if cycle[0] in adj[cycle[-1]] else None


def monocycle_numberings(cycle, elements):
    """Numberings (atom -> locant) of a mancude monocycle: heteroatoms lowest, in order of seniority."""
    n = len(cycle)
    options = []
    for start in range(n):
        for step in (1, -1):
            order = [cycle[(start + step * k) % n] for k in range(n)]
            numbering = {atom: k + 1 for k, atom in enumerate(order)}
            key = tuple(sorted(numbering[a] for a in cycle if elements[a] != "C"))
            senior = tuple(numbering[a] for el in _PRIORITY for a in cycle if elements[a] == el)
            options.append((key, senior, numbering))
    best = min((o[0], o[1]) for o in options)
    return [o[2] for o in options if (o[0], o[1]) == best]


def name_fused_monocycle(shared, comp, view_elem, view_adj, view_order):
    """Prefix for a ring made of the shared skeleton atoms and the component atoms; locants are for primes."""
    atoms = set(comp) | set(shared)
    cycle = ring_cycle(atoms, view_adj)
    if cycle is None:
        raise UnsupportedStructure("an added component that is not a simple ring")
    hetero = [view_elem[a] for a in comp if view_elem[a] != "C"]
    name, cite_hetero = _hw_name(len(cycle), [view_elem[a] for a in atoms if view_elem[a] != "C"])
    return cycle, name, cite_hetero or (len(hetero) > 1)
