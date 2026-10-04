"""Skeletal modifications of stereoparents (P-101.3): 'nor', 'homo' and 'seco' applied to a numbered graph."""

from itertools import combinations

from ._np_core import loc_key


class Skel:
    """A parent graph after modifications: atoms keep their locants, new atoms get provisional labels."""

    def __init__(self, parent, elem, adj, order, ops):
        self.parent = parent
        self.elem = elem
        self.adj = adj
        self.order = order
        self.ops = ops

    @classmethod
    def of(cls, parent):
        return cls(
            parent,
            dict(parent.elem),
            {a: set(n) for a, n in parent.adj.items()},
            dict(parent.bond_order),
            (),
        )

    def copy(self, op):
        return Skel(
            self.parent,
            dict(self.elem),
            {a: set(n) for a, n in self.adj.items()},
            dict(self.order),
            self.ops + (op,),
        )

    def remove_bond(self, a, b):
        self.adj[a].discard(b)
        self.adj[b].discard(a)
        self.order.pop(frozenset((a, b)), None)

    def add_bond(self, a, b, order=1):
        self.adj[a].add(b)
        self.adj[b].add(a)
        self.order[frozenset((a, b))] = order

    def ring_atoms(self):
        """Atoms on cycles (iteratively strip degree-1 atoms, then keep paths between cycles)."""
        degree = {a: len(n) for a, n in self.adj.items()}
        alive = set(self.adj)
        stack = [a for a in alive if degree[a] <= 1]
        while stack:
            a = stack.pop()
            if a not in alive:
                continue
            alive.discard(a)
            for n in self.adj[a]:
                if n in alive:
                    degree[n] -= 1
                    if degree[n] <= 1:
                        stack.append(n)
        return alive


def connectors(parent):
    """(atomic connectors, terminal segments, bond connectors) of a parent graph (P-101.2.5).

    Atomic connectors and terminal segments are lists of locants from the junction outward by increasing locant."""
    adj = parent.adj
    junction = {a for a, n in adj.items() if len(n) >= 3}
    atomic, terminal, seen = [], [], set()
    for start in sorted(adj, key=loc_key):
        if start in junction or start in seen or len(adj[start]) > 2:
            continue
        chain, frontier = [start], [start]
        seen.add(start)
        ends = []
        while frontier:
            a = frontier.pop()
            for n in adj[a]:
                if n in junction:
                    ends.append((a, n))
                elif n not in seen:
                    seen.add(n)
                    chain.append(n)
                    frontier.append(n)
        terminals = [a for a in chain if len(adj[a]) == 1]
        if terminals:
            terminal.append(sorted(chain, key=loc_key))
        elif ends:
            atomic.append(sorted(chain, key=loc_key))
    bond_connectors = sorted(
        (tuple(sorted(b, key=loc_key)) for b in parent.bond_order if all(a in junction for a in b)),
        key=lambda b: (loc_key(b[0]), loc_key(b[1])),
    )
    return atomic, terminal, bond_connectors


def _chain_ends(parent, chain):
    chain_set = set(chain)
    outside = [(a, n) for a in chain for n in parent.adj[a] if n not in chain_set]
    return outside


def nor_variants(parent):
    """[(op, atom)] removable atoms: the highest-locant atom of each connector, the free end of each terminal segment."""
    atomic, terminal, _ = connectors(parent)
    atoms = []
    for chain in atomic:
        atoms.append(chain[-1])
    for chain in terminal:
        end = next(a for a in chain if len(parent.adj[a]) == 1)
        atoms.append(end)
    return atoms


def apply_nor(skel, atom):
    new = skel.copy(("nor", atom))
    neighbors = sorted(new.adj[atom], key=loc_key)
    for n in neighbors:
        new.remove_bond(atom, n)
    del new.adj[atom]
    del new.elem[atom]
    if len(neighbors) == 2:
        new.add_bond(neighbors[0], neighbors[1], 1)
    return new


def homo_sites(parent):
    """[(kind, data)] insertion sites: ('atomic', (atom, outer neighbour)), ('terminal', atom), ('bond', (a, b))."""
    atomic, terminal, bond_connectors = connectors(parent)
    sites = []
    for chain in atomic:
        top = chain[-1]
        outside = [n for n in parent.adj[top] if n not in chain]
        if outside:
            sites.append(("atomic", (top, sorted(outside, key=loc_key)[-1])))
    for chain in terminal:
        end = next(a for a in chain if len(parent.adj[a]) == 1)
        sites.append(("terminal", end))
    for bond in bond_connectors:
        sites.append(("bond", bond))
    return sites


_COUNTER = [0]


def _far_end(skel, start, target):
    """The atom next to `target` on the (homo-lengthened) path from `start`."""
    current, seen = start, {start}
    while target not in skel.adj[current]:
        step = next((n for n in skel.adj[current] if n.startswith("h") and n not in seen), None)
        if step is None:
            return None
        seen.add(step)
        current = step
    return current


def apply_homo(skel, site):
    _COUNTER[0] += 1
    label = f"h{_COUNTER[0]}"
    kind, data = site
    new = skel.copy(("homo", kind, data, label))
    new.elem[label] = "C"
    new.adj[label] = set()
    if kind == "terminal":
        end = data
        while True:
            step = next((n for n in new.adj[end] if n.startswith("h") and n != label), None)
            if step is None:
                break
            end = step
        new.add_bond(end, label, 1)
        return new
    a, b = data
    tail = _far_end(new, a, b)
    order = new.order.get(frozenset((tail, b)), 1)
    new.remove_bond(tail, b)
    new.add_bond(tail, label, 1)
    new.add_bond(label, b, order if kind == "bond" else 1)
    return new


def apply_seco(skel, a, b):
    new = skel.copy(("seco", a, b))
    new.remove_bond(a, b)
    return new


def seco_bonds(parent):
    return sorted(
        (tuple(sorted(b, key=loc_key)) for b in parent.ring_bonds),
        key=lambda b: (loc_key(b[0]), loc_key(b[1])),
    )


def variants(parent, cost=1, terminal_only=False):
    """Skel graphs reachable by exactly `cost` nor/homo/seco operations."""
    max_cost = cost
    layers = [[Skel.of(parent)]]
    nors, sites, secos = nor_variants(parent), homo_sites(parent), seco_bonds(parent)
    if terminal_only:
        nors = [a for a in nors if len(parent.adj[a]) == 1]
        sites = [site for site in sites if site[0] == "terminal"]
        secos = []
    for _ in range(max_cost):
        nxt = []
        for skel in layers[-1]:
            removed = {op[1] for op in skel.ops if op[0] == "nor"}
            for atom in nors:
                if atom in skel.adj and atom not in removed and len(skel.adj[atom]) <= 2:
                    nxt.append(apply_nor(skel, atom))
            for site in sites:
                if _site_alive(skel, site):
                    nxt.append(apply_homo(skel, site))
            cut = {frozenset(op[1:]) for op in skel.ops if op[0] == "seco"}
            for a, b in secos:
                if b in skel.adj.get(a, ()) and frozenset((a, b)) not in cut:
                    nxt.append(apply_seco(skel, a, b))
        layers.append(_dedupe(nxt))
    return layers[cost]


def _site_alive(skel, site):
    kind, data = site
    if kind == "terminal":
        return data in skel.adj
    if kind == "atomic":
        top, outer = data
        return top in skel.adj and outer in skel.adj and _far_end(skel, top, outer) is not None
    a, b = data
    return a in skel.adj and b in skel.adj and _far_end(skel, a, b) is not None


def _op_key(op):
    if op[0] == "homo":
        return ("homo", op[1], str(op[2]))
    return (op[0],) + tuple(sorted(str(x) for x in op[1:]))


def _dedupe(skels):
    seen, out = set(), []
    for skel in skels:
        key = tuple(sorted(_op_key(op) for op in skel.ops))
        if key in seen:
            continue
        seen.add(key)
        out.append(skel)
    return out
