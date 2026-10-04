"""Mancude rings fused to a natural-product parent: fusion prefixes with primed locants (P-101.5.1)."""

from dataclasses import dataclass, field

from ._common import UnsupportedStructure
from ._np_core import FACES, loc_key
from ._numerals import multiplying_prefix

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
    (5, ("N", "N")): "imidazo", (6, ("N", "N")): "pyrazino",
}
_PRIMES = ["", "′", "″", "‴"]
_DONORS = ("O", "S", "Se", "Te")


@dataclass
class Fused:
    name: str
    hetero_locants: list
    attached: list
    parent_locants: list
    indicated: list
    fusion_h: list = field(default_factory=list)
    cite_attached: bool = True
    prime: int = 0
    fusion_locs: tuple = ()

    @property
    def alpha(self):
        return self.name


def hw_prefix(size, hetero):
    """Fusion prefix of a mancude monocycle and whether its heteroatom locants are cited (P-25.3.2.1)."""
    elements = tuple(sorted(hetero, key=_PRIORITY.index))
    if (size, elements) in _RETAINED:
        return _RETAINED[(size, elements)], False
    if not elements:
        if size not in _CARBO:
            raise UnsupportedStructure("a carbocyclic ring size without a fusion prefix")
        return _CARBO[size], False
    stem = _STEM.get(size)
    if stem is None:
        raise UnsupportedStructure("a ring size without a Hantzsch-Widman stem")
    counts = {}
    for el in elements:
        counts[el] = counts.get(el, 0) + 1
    multiplier = {2: "di", 3: "tri", 4: "tetra"}
    prefix = ""
    for el in sorted(counts, key=_PRIORITY.index):
        prefix += (multiplier[counts[el]] if counts[el] > 1 else "") + _HETERO_PREFIX[el]
    if size == 3 and elements == ("N",):
        text = "azirine"
    else:
        if prefix.endswith("a") and stem[0] in "aeiou":
            prefix = prefix[:-1]
        text = prefix + stem
    text = text.replace("aa", "a")
    return text[:-1] + "o" if text.endswith("e") else text + "o", len(elements) > 1


def ring_cycle(atoms, adj):
    atoms = set(atoms)
    start = sorted(atoms)[0]
    cycle, previous = [start], None
    while len(cycle) < len(atoms):
        step = next((n for n in sorted(adj[cycle[-1]]) if n in atoms and n != previous and n not in cycle), None)
        if step is None:
            return None
        previous = cycle[-1]
        cycle.append(step)
    return cycle if cycle[0] in adj[cycle[-1]] else None


def _mono_numberings(cycle, elem):
    n = len(cycle)
    options = []
    for start in range(n):
        for step in (1, -1):
            order = [cycle[(start + step * k) % n] for k in range(n)]
            numbering = {atom: k + 1 for k, atom in enumerate(order)}
            key = tuple(sorted(numbering[a] for a in cycle if elem[a] != "C"))
            senior = tuple(numbering[a] for el in _PRIORITY for a in cycle if elem[a] == el)
            options.append((key, senior, numbering))
    best = min((o[0], o[1]) for o in options)
    return [o[2] for o in options if (o[0], o[1]) == best]


def _naphthalene_numberings(atoms, adj):
    degree = {a: sum(1 for n in adj[a] if n in atoms) for a in atoms}
    fusion = [a for a in atoms if degree[a] == 3]
    if len(fusion) != 2 or len(atoms) != 10:
        return []
    rim = {a: [n for n in adj[a] if n in atoms and degree[n] == 2] for a in atoms}
    numberings = []
    for f in fusion:
        for first in [n for n in adj[f] if n in atoms and degree[n] == 2]:
            order, previous = [first], f
            while len(order) < 4:
                nxt = [n for n in adj[order[-1]] if n in atoms and n != previous and degree[n] == 2]
                if not nxt:
                    break
                previous = order[-1]
                order.append(nxt[0])
            if len(order) != 4:
                continue
            fa = next((n for n in adj[order[-1]] if n in atoms and degree[n] == 3 and n != f), None)
            if fa is None:
                continue
            tail, previous = [], fa
            current = next((n for n in adj[fa] if n in atoms and degree[n] == 2 and n not in order), None)
            while current is not None and len(tail) < 4:
                tail.append(current)
                nxt = [n for n in adj[current] if n in atoms and n != previous and degree[n] == 2]
                previous, current = current, (nxt[0] if nxt else None)
            if len(tail) != 4:
                continue
            numbering = {atom: i + 1 for i, atom in enumerate(order)}
            numbering[fa] = "4a"
            for i, atom in enumerate(tail):
                numbering[atom] = 5 + i
            numbering[f] = "8a"
            numberings.append(numbering)
    return numberings


def _locant_text(number, prime):
    return f"{number}{_PRIMES[prime]}"


def name_fused(comp, cand, view, final, parent_centers, hfaces):
    """Describe one ring fused to the skeleton at a bond (or two) of the parent."""
    image = {a: loc for loc, a in cand.mapping.items()}
    attachments = sorted({a for a, _ in comp.links})
    if len(attachments) != 2:
        raise UnsupportedStructure("an added ring attached at other than two skeleton atoms")
    sa, sb = attachments
    skeleton_adj = {
        atom: {n for n in view.adj[atom] if n in image and image[n] in cand.skel.adj[image[atom]]} for atom in image
    }
    path = _shortest(skeleton_adj, sa, sb)
    if path is None or len(path) > 3:
        raise UnsupportedStructure("an added ring that spans more than three skeleton atoms")
    ring = set(comp.atoms) | set(path)
    bonds = [(a, b) for a in ring for b in view.adj[a] if b in ring and a < b]
    rings = len(bonds) - len(ring) + 1
    if rings == 1:
        cycle = ring_cycle(ring, view.adj)
        if cycle is None:
            raise UnsupportedStructure("an added component that is not a simple ring")
        elem = {a: view.elem[a] for a in ring}
        hetero = [elem[a] for a in comp.atoms if elem[a] != "C"]
        name, cite = hw_prefix(len(cycle), hetero)
        numberings = _mono_numberings(cycle, elem)
        hetero_locants = None
    elif rings == 2:
        if any(view.elem[a] != "C" for a in ring):
            raise UnsupportedStructure("a heterocyclic bicyclic added component")
        numberings = _naphthalene_numberings(ring, view.adj)
        if not numberings:
            raise UnsupportedStructure("an added bicyclic component other than naphthalene")
        name, cite, elem = "naphtho", True, {a: "C" for a in ring}
    else:
        raise UnsupportedStructure("an added ring system of more than two rings")

    def fusion_key(numbering):
        return tuple(sorted(int(str(numbering[a]).rstrip("ab")) + (0.5 if str(numbering[a]).endswith(("a", "b")) else 0) for a in path))

    ordered_path = sorted(path, key=lambda a: loc_key(final(image[a])))

    def number_value(numbering, atom):
        text = str(numbering[atom])
        return int(text.rstrip("ab")) + (0.5 if text.endswith(("a", "b")) else 0)

    chosen = min(numberings, key=lambda n: (fusion_key(n), tuple(number_value(n, a) for a in ordered_path)))
    hetero_locants = sorted(chosen[a] for a in comp.atoms if view.elem[a] != "C") if cite and rings == 1 else []
    sp3 = _indicated(comp, path, ring, view, chosen)
    return Fused(
        name=name,
        hetero_locants=hetero_locants,
        attached=[chosen[a] for a in ordered_path],
        parent_locants=[final(image[a]) for a in ordered_path],
        indicated=sp3[0],
        fusion_h=[(final(image[a]), a) for a in ordered_path if _fusion_h(a, view, ring, image, parent_centers)],
        cite_attached=name not in _CARBO.values() or rings == 2,
        fusion_locs=tuple(image[a] for a in path),
    )


def _shortest(adj, start, goal):
    frontier, parents = [start], {start: None}
    while frontier:
        nxt = []
        for a in frontier:
            for n in adj[a]:
                if n not in parents:
                    parents[n] = a
                    nxt.append(n)
        frontier = nxt
    if goal not in parents:
        return None
    path = [goal]
    while parents[path[-1]] is not None:
        path.append(parents[path[-1]])
    return path[::-1]


def _double_in(atom, ring, view):
    return any(view.order.get(frozenset((atom, n)), 1) == 2 for n in view.adj[atom] if n in ring)


def _indicated(comp, path, ring, view, numbering):
    """Component locants carrying indicated hydrogen (P-101.5.1.2); raises when hydro prefixes would be needed."""
    donors = {a for a in comp.atoms if view.elem[a] in _DONORS or (view.elem[a] == "N" and not _double_in(a, ring, view))}
    nondonor = [a for a in comp.atoms if a not in donors]
    sp2_fusion = sum(1 for a in path if _double_in(a, ring, view))
    need = (len(nondonor) + sp2_fusion) % 2
    saturated = [a for a in nondonor if not _double_in(a, ring, view)]
    if len(saturated) != need:
        raise UnsupportedStructure("a fused component needing hydro prefixes")
    return [numbering[a] for a in saturated], donors


def _fusion_h(atom, view, ring, image, parent_centers):
    if _double_in(atom, ring, view) or image[atom] in parent_centers:
        return False
    return view.mol.GetAtomWithIdx(atom).GetTotalNumHs() > 0


def assign_primes(fused):
    """Order the fused components for citation and give each its prime (P-101.5.1.1)."""
    ordered = []
    groups = {}
    for f in fused:
        groups.setdefault((f.name, tuple(f.hetero_locants)), []).append(f)
    for key in sorted(groups, key=lambda k: k[0]):
        ordered.append(groups[key])
    prime = 0
    for items in ordered:
        for f in items:
            prime += 1
            f.prime = prime
    return ordered


def fusion_text(groups):
    text = ""
    for items in groups:
        brackets = []
        for f in items:
            attached = ",".join(_locant_text(n, f.prime) for n in f.attached) if f.cite_attached else ""
            parent = ",".join(f.parent_locants)
            brackets.append(f"{attached}:{parent}" if attached else parent)
        hetero_text = f"[{','.join(str(h) for h in items[0].hetero_locants)}]" if items[0].hetero_locants else ""
        multiplier = "" if len(items) == 1 else multiplying_prefix(len(items), compound=bool(hetero_text))
        text += f"{multiplier}{hetero_text}{items[0].name}[{';'.join(brackets)}]"
    return text


def indicated_texts(groups):
    found = []
    for items in groups:
        for f in items:
            found += [(f.prime, n) for n in f.indicated]
    return [f"{n}{_PRIMES[prime]}H" for prime, n in sorted(found)]
