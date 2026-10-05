"""P-25.4 bridged fused ring systems: a fused ring system (P-25.1-25.3) plus independent and dependent bridges (divalent,
polyvalent, composite, with ring and fused-ring units), named with the bridge prefixes cited in front of the fused parent,
the attachment locants of the bridges, and the bridge atoms numbered on from the fused system (P-25.4.3-5)."""

from dataclasses import dataclass, field
from itertools import groupby, product

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure
from ._fused_numbering import HETERO_RANK, FusedSystem, _locant_key
from ._fusion_components import identify, skeleton
from ._fusion_name import Context, _fusion_name_core, system_numbering_options
from ._numerals import alkane_name, numerical_term

_MAX_CYCLE = 12
_MAX_STATES = 6000
_MAX_CLUSTER = 14
_MAX_PARTS = 4
_MAX_DECOMPOSITIONS = 200
_PHANE_NODES = 6
_NO_KEY = (-1, "", 0)
_HETERO_WORDS = {"S": "sulfano", "Se": "selano", "Te": "tellano", "Si": "silano", "P": "phosphano", "B": "borano", "Sn": "stannano"}
_HYDRIDES = {
    "N": "azane", "O": "oxidane", "S": "sulfane", "Se": "selane", "Te": "tellane", "P": "phosphane", "As": "arsane",
    "Sb": "stibane", "Si": "silane", "B": "borane", "Ge": "germane", "Sn": "stannane",
}
_POLYVALENT = {
    ("C", (1, 2)): "metheno", ("C", (1, 1, 1)): "epimethanetriyl", ("N", (1, 2)): "azeno", ("N", (1, 1, 1)): "epinitrilo",
    ("P", (1, 2)): "phospheno", ("P", (1, 1, 1)): "epiphosphanetriyl", ("Si", (1, 1, 1, 1)): "episilanetetrayl",
}


@dataclass
class BridgedParent:
    stem: str
    position_of: dict
    capable: frozenset
    consumed: frozenset
    bridge_atoms: frozenset
    key: tuple


@dataclass
class Token:
    """A free valence of a bridge unit: the bond from `inner` (in the unit) to `outer`; kind L (inside the bridge), F, B."""

    inner: int
    outer: int
    order: int
    kind: str


@dataclass
class Part:
    atoms: frozenset
    units: list
    independent: bool
    valence: int
    fused_bonds: list = field(default_factory=list)


def _fused_valid(mol, graph, cycles, subset):
    atoms, edges = set(), set()
    for i in subset:
        c = cycles[i]
        atoms |= set(c)
        edges |= {frozenset((c[k], c[(k + 1) % len(c)])) for k in range(len(c))}
    for a, b in graph.subgraph(atoms).edges():
        if frozenset((a, b)) not in edges:
            return None
    order = sorted(atoms)
    index = {a: i for i, a in enumerate(order)}
    sk = skeleton([mol.GetAtomWithIdx(a).GetSymbol() for a in order], [tuple(index[x] for x in e) for e in edges])
    if sk.GetRingInfo().NumRings() != len(subset):
        return None
    try:
        FusedSystem(sk)
    except UnsupportedStructure:
        return None
    return order, sk


def _candidates(mol, atoms):
    graph = nx.Graph()
    graph.add_nodes_from(atoms)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in atoms and b in atoms:
            graph.add_edge(a, b)
    cycles = [tuple(c) for c in nx.simple_cycles(graph, length_bound=_MAX_CYCLE)]
    edge_sets = [{frozenset((c[k], c[(k + 1) % len(c)])) for k in range(len(c))} for c in cycles]
    total = graph.number_of_edges() - graph.number_of_nodes() + 1
    seen, valid, by_edges = set(), {}, {}
    frontier = {frozenset([i]) for i in range(len(cycles))}
    while frontier and len(seen) < _MAX_STATES:
        following = set()
        for state in frontier:
            if state in seen:
                continue
            seen.add(state)
            result = _fused_valid(mol, graph, cycles, state)
            if result is None:
                continue
            used = set().union(*[edge_sets[i] for i in state])
            size = sum(len(cycles[i]) for i in state)
            key = frozenset(used)
            if key in by_edges and by_edges[key][1] <= size:
                continue
            if key in by_edges:
                del valid[by_edges[key][0]]
            by_edges[key] = (state, size)
            valid[state] = result
            if len(state) >= total:
                continue
            for j in range(len(cycles)):
                if j not in state and edge_sets[j] & used:
                    following.add(state | {j})
        frontier = following
    return graph, valid


def _order(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _chain_word(orders):
    """Word and double-bond locants of a hydrocarbon bridge read in the given direction."""
    n = len(orders) + 1
    if not any(o != 1.0 for o in orders):
        return alkane_name(n)[:-1] + "o", ()
    if any(o not in (1.0, 2.0) for o in orders):
        raise UnsupportedStructure("a triple bond in a hydrocarbon bridge")
    locants = tuple(i + 1 for i, o in enumerate(orders) if o == 2.0)
    if n == 2:
        return "etheno", locants
    stem = alkane_name(n)[:-3] + ("a" if len(locants) > 1 else "")
    word = "eno" if len(locants) == 1 else numerical_term(len(locants)) + "eno"
    return f"{stem}[{','.join(map(str, locants))}]{word}", locants


_STANDARD_BONDING = {"N": 3, "P": 3, "As": 3, "Sb": 3, "Bi": 3, "B": 3, "O": 2, "S": 2, "Se": 2, "Te": 2, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4}


def _hetero_word(mol, atoms, first):
    symbol = mol.GetAtomWithIdx(atoms[0]).GetSymbol()
    count = len(atoms)
    multiplier = numerical_term(count) if count > 1 else ""
    odd = {mol.GetAtomWithIdx(a).GetTotalValence() for a in atoms if mol.GetAtomWithIdx(a).GetTotalValence() != _STANDARD_BONDING.get(symbol)}
    if odd:
        if count > 1 or symbol not in _HETERO_WORDS:
            raise UnsupportedStructure("a heteroatom of nonstandard bonding number in a bridge chain")
        return f"λ{odd.pop()}-" + _HETERO_WORDS[symbol]  # P-25.4.2.1.4; the hyphen is not shown in an example
    if symbol == "N":
        orders = [_order(mol, a, b) for a, b in zip(atoms, atoms[1:])]
        if count == 1:
            return "azano"
        if count == 2:
            return "diazeno" if orders[0] == 2.0 else "diazano"
        raise UnsupportedStructure("a nitrogen chain in a bridge")
    if symbol == "O":
        return ("ep" if count == 1 else "epi") * first + multiplier + "oxy"
    if symbol in _HETERO_WORDS:
        return multiplier + _HETERO_WORDS[symbol]
    raise UnsupportedStructure(f"a {symbol} bridge")


_ring_cache = {}


def _ring_info(mol, ring_atoms):
    """Name, fusion prefix and numberings (atom -> locant) of the ring or fused ring system a bridge unit is."""
    order = sorted(ring_atoms)
    index = {a: i for i, a in enumerate(order)}
    edges = [
        (index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()])
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in index and b.GetEndAtomIdx() in index
    ]
    symbols = [mol.GetAtomWithIdx(a).GetSymbol() for a in order]
    key = (tuple(symbols), tuple(sorted(edges)))
    if key not in _ring_cache:
        sk = skeleton(symbols, edges)
        if sk.GetRingInfo().NumRings() == 1:
            found = [c for c in identify(sk) if len(c.ring_sizes) == 1]
            if not found:
                raise UnsupportedStructure("an unnamed cyclic bridge")
            _ring_cache[key] = (found[0].name, found[0].prefix, found[0].numberings)
        else:
            name, root = _fusion_name_core(sk)
            known = [c.prefix for c in identify(sk) if c.name == name]
            _ring_cache[key] = (name, known[0] if known else "", system_numbering_options(Context(sk), name, root))
    name, prefix, numberings = _ring_cache[key]
    return name, prefix, [{order[i]: loc for i, loc in n.items()} for n in numberings]


def _ring_word(name, prefix, first):
    word = name[:-1] + "o" if name.endswith("e") else name + "o"
    if name == "benzene":
        return "benzeno"
    if word == prefix:
        return ("ep" if word[0] in "io" else "epi") + word if first else word
    return word


def _units_of(mol, graph, comp):
    """The path of units (acyclic runs of one element, rings, fused ring systems) a bridge is made of."""
    sub = graph.subgraph(comp)
    groups = []
    for cycle in nx.cycle_basis(sub):
        merged = set(cycle)
        for group in groups[:]:
            if group & merged:
                merged |= group
                groups.remove(group)
        groups.append(merged)
    node_of = {a: ("ring", frozenset(g)) for g in groups for a in g}
    for a in comp:
        node_of.setdefault(a, ("atom", a))
    quotient = nx.Graph()
    quotient.add_nodes_from(set(node_of.values()))
    for a, b in sub.edges():
        if node_of[a] != node_of[b]:
            quotient.add_edge(node_of[a], node_of[b])
    if not nx.is_connected(quotient) or any(d > 2 for _, d in quotient.degree()) or quotient.number_of_edges() != quotient.number_of_nodes() - 1:
        raise UnsupportedStructure("a branched or cyclic bridge")
    ends = [n for n in quotient if quotient.degree(n) <= 1]
    path = nx.shortest_path(quotient, ends[0], ends[-1]) if len(ends) > 1 else ends
    units = []
    for kind, payload in path:
        if kind == "ring":
            _ring_info(mol, payload)
            units.append({"kind": "ring", "atoms": sorted(payload)})
            continue
        symbol = mol.GetAtomWithIdx(payload).GetSymbol()
        if units and units[-1]["kind"] != "ring" and units[-1]["symbol"] == symbol:
            units[-1]["atoms"].append(payload)
        else:
            units.append({"kind": "chain" if symbol == "C" else "hetero", "atoms": [payload], "symbol": symbol})
    ring_units = sum(1 for u in units if u["kind"] == "ring")
    chain_atoms = sum(len(u["atoms"]) for u in units if u["kind"] != "ring")
    if ring_units and chain_atoms + ring_units + 1 >= _PHANE_NODES:
        raise UnsupportedStructure("a bridge that makes a cyclophane (P-52.2.5.1)")
    return units


def _connected_subsets(sub):
    """All connected vertex subsets of a small graph."""
    found, frontier = set(), {frozenset([v]) for v in sub}
    while frontier:
        found |= frontier
        frontier = {s | {n} for s in frontier for v in s for n in sub[v] if n not in s} - found
    return found


def _decompositions(mol, graph, fused, cluster, cache, fused_rings):
    """Ways to read a connected group of non-fused atoms as independent and dependent bridges (P-25.4.1.8-9)."""
    if len(cluster) > _MAX_CLUSTER:
        raise UnsupportedStructure("a bridge system with too many atoms")
    sub = graph.subgraph(cluster)
    key = ("partitions", frozenset(cluster))
    if key not in cache:
        units = {}
        for subset in _connected_subsets(sub):
            if subset not in cache:
                try:
                    cache[subset] = _units_of(mol, graph, subset)
                except UnsupportedStructure:
                    cache[subset] = None
            if cache[subset] is not None:
                units[subset] = cache[subset]
        partitions = []

        def cover(left, chosen):
            if not left:
                partitions.append(list(chosen))
                return
            if len(chosen) == _MAX_PARTS:
                return
            first = min(left)
            for subset in units:
                if first in subset and subset <= left:
                    cover(left - subset, chosen + [subset])

        cover(frozenset(cluster), [])
        whole_rings = [set(r) for r in cache["aromatic_rings"] if set(r) <= cluster]
        partitions = [p for p in partitions if all(any(r <= part for part in p) for r in whole_rings)]
        cache[key] = partitions
    partitions = cache[key]
    found = []
    for comps in partitions:
        owner = {a: i for i, c in enumerate(comps) for a in c}
        fused_bonds = [[(b, n, _order(mol, b, n)) for b in c for n in graph[b] if n in fused] for c in comps]
        links = [[] for _ in comps]
        for a, b in sub.edges():
            if owner[a] != owner[b]:
                links[owner[a]].append((owner[b], _order(mol, a, b)))
                links[owner[b]].append((owner[a], _order(mol, a, b)))
        for roles in product((True, False), repeat=len(comps)):
            if not _roles_valid(roles, fused_bonds, links):
                continue
            if any(r and _forms_fused_ring(graph, fused, fb, fused_rings) for r, fb in zip(roles, fused_bonds)):
                continue
            parts = []
            for i, c in enumerate(comps):
                valence = sum(o for _, _, o in fused_bonds[i]) + (0 if roles[i] else sum(o for _, o in links[i]))
                parts.append(Part(c, cache[c], roles[i], int(valence), fused_bonds[i]))
            found.append(parts)
    if not found:
        raise UnsupportedStructure("this bridge system has no bridge reading")
    return found


def _forms_fused_ring(graph, fused, fused_bonds, fused_rings):
    """P-25.4.1.2 (b): a divalent bridge that closes a new ortho- or peri-fused ring is a fusion component, not a bridge;
    only the bond common to two rings (d) may be bridged."""
    ends = [f for _, f, _ in fused_bonds]
    if len(ends) != 2 or ends[0] == ends[1]:
        return False
    a, b = ends
    both = sum(1 for r in fused_rings if a in r and b in r)
    if graph.has_edge(a, b):
        return both < 2
    return any(
        sum(c in r for r in fused_rings) >= 2 and not both
        for c in (set(graph[a]) & set(graph[b]) & fused)
    )


def _roles_valid(roles, fused_bonds, links):
    for i, independent in enumerate(roles):
        if independent:
            if len(fused_bonds[i]) < 2 or any(roles[j] for j, _ in links[i]):
                return False
        elif not links[i] or any(not roles[j] for j, _ in links[i]) or len(fused_bonds[i]) + len(links[i]) < 2:
            return False
    return True


def _metrics(parts):
    """P-25.4.3.4.2 (e)-(h): fewer polyvalent bridges, fewer dependent bridges, fewer atoms in dependent bridges, more divalent bridges.
    The P-25.4.1.9 example (epimethanetriyl + epipropane[1,2,3]triyl) contradicts (e) and P-25.4.3.3 (a); the criteria win."""
    dependent = [p for p in parts if not p.independent]
    return (
        sum(1 for p in parts if p.valence >= 3),
        len(dependent),
        sum(len(p.atoms) for p in dependent),
        -sum(1 for p in parts if p.valence == 2),
    )


def _tokens(mol, atoms, fused, part, unit_atoms):
    unit = set(unit_atoms)
    out = []
    for x in unit_atoms:
        for n in mol.GetAtomWithIdx(x).GetNeighbors():
            y = n.GetIdx()
            if y not in atoms or y in unit:
                continue
            kind = "L" if y in part.atoms else "F" if y in fused else "B"
            if kind == "B" and part.independent:
                continue
            order = _order(mol, x, y)
            if order not in (1.0, 2.0):
                raise UnsupportedStructure("a triple bond or aromatic bond at a bridge")
            out.append(Token(x, y, int(order), kind))
    return out


def _multiplied(n):
    return "" if n == 1 else numerical_term(n)


def _ending(singles, doubles, show):
    pieces = []
    for locants, word in ((singles, "yl"), (doubles, "ylidene")):
        if locants:
            pieces.append(((f"[{','.join(map(str, locants))}]" if show else ""), _multiplied(len(locants)) + word))
    return "".join(a + b for a, b in pieces), pieces[0][1][0]


def _substituent_word(mol, seq, singles, doubles, first):
    """P-25.4.2.2.2: a polyvalent polyatomic bridge named as the polyvalent substituent group."""
    n = len(seq)
    symbol = mol.GetAtomWithIdx(seq[0]).GetSymbol()
    inner = [_order(mol, a, b) for a, b in zip(seq, seq[1:])]
    if symbol == "C":
        if any(o not in (1.0, 2.0) for o in inner):
            raise UnsupportedStructure("a triple bond in a polyvalent bridge")
        ene = [i + 1 for i, o in enumerate(inner) if o == 2.0]
        if n == 1:
            stem = "methane"
        elif not ene:
            stem = alkane_name(n)
        elif n == 2:
            stem = "ethene"
        else:
            stem = alkane_name(n)[:-3] + ("a" if len(ene) > 1 else "") + f"[{','.join(map(str, ene))}]" + ("ene" if len(ene) == 1 else numerical_term(len(ene)) + "ene")
    else:
        if symbol not in _HYDRIDES or any(o != 1.0 for o in inner):
            raise UnsupportedStructure(f"a polyvalent {symbol} bridge")
        ene = []
        stem = (numerical_term(n) if n > 1 else "") + _HYDRIDES[symbol]
    show = not (n == 1 or (n == 2 and sorted(singles + doubles) == [1, 2]))
    text, letter = _ending(singles, doubles, show)
    if letter == "y":
        stem = stem[:-1]
    return ("epi" if first else "") + stem + text, tuple(ene)


def _run_word(mol, seq, tokens, first, partner):
    """Word of an acyclic run (carbon or one heteroelement), the tokens in the order the name cites them, and its rank."""
    n = len(seq)
    pos = {a: i + 1 for i, a in enumerate(seq)}
    symbol = mol.GetAtomWithIdx(seq[0]).GetSymbol()

    def listed(order):
        return sorted(((pos[t.inner], t) for t in tokens if t.order == order), key=lambda e: (e[0], partner(e[1].outer) if e[1].kind != "L" else _NO_KEY))

    singles, doubles = listed(1), listed(2)
    where = sorted(p for p, _ in singles)
    shape = tuple(sorted(t.order for t in tokens))
    if not doubles and len(singles) == 2 and where == ([1, n] if n > 1 else [1, 1]):
        inner = [_order(mol, a, b) for a, b in zip(seq, seq[1:])]
        word, ene = _chain_word(inner) if symbol == "C" else (_hetero_word(mol, seq, first), ())
        locs = ((1, n), (), ene)
    elif n == 1 and (symbol, shape) in _POLYVALENT:
        word = _POLYVALENT[(symbol, shape)]
        if not first and word.startswith("epi"):
            word = word[3:]
        locs = ((1,) * len(singles), (1,) * len(doubles), ())
    else:
        word, ene = _substituent_word(mol, seq, [p for p, _ in singles], [p for p, _ in doubles], first)
        locs = (tuple(p for p, _ in singles), tuple(p for p, _ in doubles), ene)
    cited = [t for _, t in singles + doubles if t.kind != "L"]
    rank = (2, 0, -n, word, locs) if symbol == "C" else (0, HETERO_RANK.get(symbol, 99), 0, word, locs)
    return word, cited, rank


def _ring_unit_word(mol, unit, tokens, first, partner, side):
    name, prefix, numberings = _ring_info(mol, unit["atoms"])
    ordered = sorted(tokens, key=lambda t: (side(t), t.inner))
    singles = [t for t in ordered if t.order == 1]
    doubles = [t for t in ordered if t.order == 2]
    divalent = not doubles and len(singles) == 2
    best = None
    for numbering in numberings:
        if divalent:
            locs = [numbering[t.inner] for t in ordered]
            key = (tuple(sorted(_locant_key(str(l)) for l in locs)), tuple(_locant_key(str(l)) for l in locs))
            groups = [locs]
        else:
            groups = [[numbering[t.inner] for t in sorted(g, key=lambda t: _locant_key(str(numbering[t.inner])))] for g in (singles, doubles)]
            key = tuple(tuple(_locant_key(str(l)) for l in g) for g in groups)
        if best is None or key < best[0]:
            best = (key, groups, numbering)
    key, groups, numbering = best
    if divalent:
        word = f"[{','.join(map(str, groups[0]))}]" + _ring_word(name, prefix, first)
        cited = ordered
    else:
        text, _ = _ending([str(l) for l in groups[0]], [str(l) for l in groups[1]], True)
        word = ("epi" if first else "") + name + text
        cited = sorted(ordered, key=lambda t: (t.order, _locant_key(str(numbering[t.inner])), partner(t.outer) if t.kind != "L" else _NO_KEY))
    return word, [t for t in cited if t.kind != "L"], (1, 0, 0, word, key)


def _part_options(mol, atoms, fused, part, partner):
    """The directions of the bridge that give the senior name: words, cited tokens and the unit order."""
    candidates = []
    for forward in (True, False):
        units = part.units if forward else part.units[::-1]
        words, ranks, cited = [], [], []
        for index, unit in enumerate(units):
            tokens = _tokens(mol, atoms, fused, part, unit["atoms"])
            if unit["kind"] == "ring":
                before = set(units[index - 1]["atoms"]) if index else set()
                last = len(units) - 1

                def side(t, before=before, index=index, last=last, forward=forward, single=len(units) == 1):
                    if t.kind == "L":
                        return 0 if t.outer in before else 2
                    if single:
                        return t.inner if forward else -t.inner
                    return 0 if index == 0 else 2 if index == last else 1

                word, c, rank = _ring_unit_word(mol, unit, tokens, index == 0, partner, side)
            else:
                seq = unit["atoms"] if forward else unit["atoms"][::-1]
                word, c, rank = _run_word(mol, seq, tokens, index == 0, partner)
            words.append(word)
            ranks.append(rank)
            cited.extend(c)
        candidates.append({"units": units, "words": words, "ranks": ranks, "cited": cited})
    best = min(c["ranks"] for c in candidates)
    return [c for c in candidates if c["ranks"] == best]


def _periphery(mol, atoms):
    """Atoms of a fused ring unit in peripheral order, and its fusion atoms."""
    order = sorted(atoms)
    index = {a: i for i, a in enumerate(order)}
    sk = skeleton(
        [mol.GetAtomWithIdx(a).GetSymbol() for a in order],
        [(index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()]) for b in mol.GetBonds() if b.GetBeginAtomIdx() in index and b.GetEndAtomIdx() in index],
    )
    rings = [list(r) for r in sk.GetRingInfo().AtomRings()]
    if any(sum(i in r for r in rings) > 2 for i in range(len(order))):
        raise UnsupportedStructure("a bridge ring system with interior atoms")
    count = {}
    for ring in rings:
        for k in range(len(ring)):
            e = frozenset((ring[k], ring[(k + 1) % len(ring)]))
            count[e] = count.get(e, 0) + 1
    outer = nx.Graph([tuple(e) for e, c in count.items() if c == 1])
    fusion = {order[i] for i in range(len(order)) if sum(i in r for r in rings) > 1}
    return [order[i] for i in nx.cycle_basis(outer)[0]], fusion


def _walks(mol, unit_atoms, start):
    """The orders in which the atoms of a ring unit are numbered when numbering starts at `start`."""
    if len(unit_atoms) == 1:
        return [([start], set())]
    cycle, fusion = _periphery(mol, unit_atoms)
    if start not in cycle:
        raise UnsupportedStructure("a bridge ring numbered from an interior atom")
    k = cycle.index(start)
    forward = cycle[k:] + cycle[:k]
    return [(forward, fusion), ([forward[0]] + forward[1:][::-1], fusion)]


def _bridge_numberings(mol, atoms, fused, part, position, counter):
    """P-25.4.4: bridge atoms numbered on from `counter`, beginning at the end attached to the bridgehead with the highest
    locant; every unit is numbered completely before the next; fusion atoms of ring units take letter locants."""
    candidates = []
    for forward in (True, False):
        units = part.units if forward else part.units[::-1]
        options = []
        for index, unit in enumerate(units):
            if unit["kind"] != "ring":
                options.append([(unit["atoms"] if forward else unit["atoms"][::-1], set())])
                continue
            tokens = _tokens(mol, atoms, fused, part, unit["atoms"])
            if index:
                before = set(units[index - 1]["atoms"])
                starts = sorted({t.inner for t in tokens if t.kind == "L" and t.outer in before})
            else:
                starts = sorted({t.inner for t in tokens if t.kind != "L"})
            options.append([w for s in starts for w in _walks(mol, unit["atoms"], s)])
        for combo in product(*options):
            candidates.append(([a for seq, _ in combo for a in seq], set().union(*[f for _, f in combo])))
    if not candidates:
        raise UnsupportedStructure("a bridge that cannot be numbered")

    def outside(a):
        return [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in position and n.GetIdx() not in part.atoms]

    def start_key(c):
        keys = [_locant_key(str(position[n])) for n in outside(c[0][0])]
        return max(keys) if keys else _NO_KEY

    top = max(start_key(c) for c in candidates)
    candidates = [c for c in candidates if start_key(c) == top]

    def rest_key(c):
        sequence = c[0]
        return (
            [i for i, a in enumerate(sequence) if mol.GetAtomWithIdx(a).GetSymbol() != "C"],
            [i for i, a in enumerate(sequence) if outside(a)],
        )

    best = min(rest_key(c) for c in candidates)
    out, seen = [], set()
    for sequence, fusion in candidates:
        if rest_key((sequence, fusion)) != best:
            continue
        labels, number, run = {}, counter, 0
        for a in sequence:
            if a in fusion and a != sequence[0]:
                run += 1
                labels[a] = f"{number}{chr(ord('a') + run - 1)}"
            else:
                number, run = number + 1, 0
                labels[a] = number
        signature = tuple(sorted(labels.items()))
        if signature not in seen:
            seen.add(signature)
            out.append((labels, number))
    return out


def _higher_number(numbering):
    best = 0
    for loc in numbering.values():
        best = max(best, _locant_key(str(loc))[0])
    return best


def _fused_parent(order, sk):
    name, root = _fusion_name_core(sk)
    return name, system_numbering_options(Context(sk), name, root), root


def _kekule(mol, single_bonds=()):
    """A Kekule form in which the aromatic bonds (begin, end) of `single_bonds` are single, if there is one."""
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Exception:
        return kekule
    single_bonds = [(a, b) for a, b in single_bonds if mol.GetBondBetweenAtoms(a, b).GetIsAromatic()]
    if not single_bonds or all(_order(kekule, a, b) == 1.0 for a, b in single_bonds):
        return kekule
    forced = Chem.RWMol(mol)
    for a, b in single_bonds:
        bond = forced.GetBondBetweenAtoms(a, b)
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    try:
        Chem.Kekulize(forced, clearAromaticFlags=True)
    except Exception:
        return None
    return forced.GetMol()


def bridged_parents(mol, atoms):
    """Names of the ring system `atoms` as a fused ring system with bridges: list of `BridgedParent`, best choice only."""
    if any(mol.GetAtomWithIdx(a).GetFormalCharge() or mol.GetAtomWithIdx(a).GetNumRadicalElectrons() for a in atoms):
        raise UnsupportedStructure("a charged or radical ring atom in a bridged fused system")
    aromatic = Chem.Mol(mol)
    mol = _kekule(aromatic)
    graph, valid = _candidates(mol, atoms)
    ranked = []
    cache = {"aromatic_rings": [r for r in aromatic.GetRingInfo().AtomRings() if all(aromatic.GetAtomWithIdx(a).GetIsAromatic() for a in r)]}

    def hetero_count(order):
        return sum(1 for a in order if mol.GetAtomWithIdx(a).GetSymbol() != "C")

    states = sorted(
        (item for item in valid.items() if len(item[0]) >= 2 and len(item[1][0]) != len(atoms) and sum(len(c) >= 5 for c in item[1][1].GetRingInfo().AtomRings()) >= 2),
        key=lambda item: (-len(item[0]), -len(item[1][0]), hetero_count(item[1][0])),
    )
    for _, group in groupby(states, key=lambda item: (-len(item[0]), -len(item[1][0]), hetero_count(item[1][0]))):
        for state, (order, sk) in group:
            fused = set(order)
            fused_rings = [{order[i] for i in ring} for ring in sk.GetRingInfo().AtomRings()]
            clusters = list(nx.connected_components(graph.subgraph(set(graph.nodes) - fused)))
            kekule = _kekule(aromatic, [(b, f) for c in clusters for b in c for f in graph[b] if f in fused]) or mol
            try:
                readings = [_decompositions(kekule, graph, fused, c, cache, fused_rings) for c in clusters]
                name, numberings, root = _fused_parent(order, sk)
            except UnsupportedStructure:
                continue
            senior = root.comp.senior_key if hasattr(root, "comp") else ()
            for combo in list(product(*readings))[:_MAX_DECOMPOSITIONS]:
                parts = sum(combo, [])
                ranked.append(((-len(state), -len(order), hetero_count(order), senior) + _metrics(parts), order, name, numberings, parts, kekule))
        if ranked:
            break
    if not ranked:
        raise UnsupportedStructure("this ring system has no fused ring system with bridges")
    ranked.sort(key=lambda r: r[0])
    out = []
    for rank, order, name, numberings, parts, kekule in ranked:
        if rank != ranked[0][0]:
            break
        try:
            out.extend(_name(kekule, atoms, order, name, numberings, parts))
        except UnsupportedStructure:
            continue
    if not out:
        raise UnsupportedStructure("no bridged name")
    best = min(p.key for p in out)
    return [p for p in out if p.key == best]


# Bridgehead CH atoms count as saturated (hydro prefixes), as in the P-25.4.3.4.2 (j) example; the PIN lines of
# P-25.4.3.2.1 (9,10-ethanoanthracene, 1,4-epoxynaphthalene) omit them and contradict that example.
def _unsaturation_key(mol, capable, consumed):
    """P-25.4.3.4.2 (j): fewer saturated positions means more noncumulative double bonds in the parent ring system."""
    live = set(capable) - set(consumed)
    return sum(
        1 for a in live if not any(mol.GetBondBetweenAtoms(a, n.GetIdx()).GetBondTypeAsDouble() == 2.0 and n.GetIdx() in live for n in mol.GetAtomWithIdx(a).GetNeighbors())
    )


def _consumed(mol, atoms, order, parts):
    """Ring atoms that give up their ring double bond: quaternary, hetero or double-bonded through a bridge bond."""
    home = {a: set(order) for a in order}
    for p in parts:
        for u in p.units:
            if u["kind"] == "ring":
                home.update({a: set(u["atoms"]) for a in u["atoms"]})
    out = set()
    for a, own in home.items():
        atom = mol.GetAtomWithIdx(a)
        neighbours = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in atoms]
        outside = [n for n in neighbours if n not in own]
        if outside and (len(neighbours) >= 4 or atom.GetSymbol() != "C" or any(_order(mol, a, n) > 1.0 for n in outside)):
            out.add(a)
    return home, frozenset(out)


def _name(mol, atoms, order, fused_name, numberings, parts):
    results = []
    for local in numberings:
        numbering = {order[i]: loc for i, loc in local.items()}
        results.extend(_assemble(mol, atoms, set(order), order, fused_name, numbering, parts))
    return results


def _number_group(mol, atoms, fused, group, position, counter):
    """Number the bridge atoms of `group` one bridge after the other, higher bridgehead first (P-25.4.5.2)."""
    def heights(part):
        tokens = [t for u in part.units for t in _tokens(mol, atoms, fused, part, u["atoms"]) if t.kind != "L"]
        return sorted((_locant_key(str(position[t.outer])) for t in tokens), reverse=True)

    def provisional(part):
        partner = lambda x: _locant_key(str(position[x])) if x in position else _NO_KEY
        return "".join(_part_options(mol, atoms, fused, part, partner)[0]["words"])

    states = [(position, counter)]
    for part in sorted(sorted(group, key=provisional), key=heights, reverse=True):
        following = []
        for current, number in states:
            for labels, last in _bridge_numberings(mol, atoms, fused, part, current, number):
                following.append(({**current, **labels}, last))
        states = following
    return states


def _assemble(mol, atoms, fused, order, fused_name, numbering, parts):
    independent = [p for p in parts if p.independent]
    dependent = [p for p in parts if not p.independent]
    results = []
    for position, counter in _number_group(mol, atoms, fused, independent, dict(numbering), _higher_number(numbering)):
        for settled, _ in _number_group(mol, atoms, fused, dependent, position, counter):
            def partner(x, settled=settled):
                return _locant_key(str(settled[x])) if x in settled else _NO_KEY

            per_part = [_part_options(mol, atoms, fused, p, partner) for p in independent + dependent]
            for chosen in product(*per_part):
                entries = []
                for part, option in zip(independent + dependent, chosen):
                    entries.append(("".join(option["words"]), [settled[t.outer] for t in option["cited"]], option, part))
                # P-25.4.3.1 alphabetical order; the P-25.4.5.1 example cites diepoxy before benzeno, against that rule
                entries.sort(key=lambda e: (e[3].independent, _citation_name(e[0]), tuple(_locant_key(str(l)) for l in e[1])))
                key = (
                    tuple(sorted(_locant_key(str(l)) for e in entries for l in e[1])),
                    tuple(tuple(_locant_key(str(l)) for l in e[1]) for e in entries),
                )
                home, consumed = _consumed(mol, atoms, order, parts)
                capable = frozenset(home)
                results.append(
                    BridgedParent(
                        _prefix_text(entries) + fused_name,
                        dict(settled),
                        capable,
                        consumed,
                        frozenset(a for p in parts for a in p.atoms),
                        key + (_unsaturation_key(mol, capable, consumed),),
                    )
                )
    return results


def _citation_name(text):
    return text.replace("[", "").lstrip("0123456789,]")


def _prefix_text(entries):
    groups = {}
    for text, locants, option, part in entries:
        groups.setdefault((not part.independent, text), []).append((locants, option, part))
    pieces = []
    for dependent, text in sorted(groups, key=lambda k: (not k[0], _citation_name(k[1]))):
        members = groups[(dependent, text)]
        located = ":".join(",".join(str(l) for l in m[0]) for m in sorted(members, key=lambda m: tuple(_locant_key(str(l)) for l in m[0])))
        composite = len(members[0][1]["units"]) > 1 or members[0][2].valence >= 3
        if len(members) > 1:
            count = len(members)
            word = ({2: "bis", 3: "tris", 4: "tetrakis"}[count] + f"({text})") if composite else numerical_term(count) + text
        else:
            word = f"({text})" if composite else text
        pieces.append(f"{located}-{word}")
    return "-".join(pieces)
