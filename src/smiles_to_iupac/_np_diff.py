"""nor/homo/seco operations read from the difference between a parent and a molecule (P-101.3.7).

The parent is cut into junctions (degree >= 3) and connectors (chains between junctions, or ending freely). The fewest
operations is a subgraph edit distance on connector lengths, searched VF2-style with arc-consistent domains. The search
starts from connected pieces that no operation can have changed, chosen through obstructions the molecule cannot contain.
The caller deepens the cost (P-101.3.7.1); the embedding that follows verifies each proposal.
"""

from rdkit import Chem

from ._np_core import loc_key
from ._np_skel import Skel, _site_alive, apply_homo, apply_nor, apply_seco, homo_sites, is_branch_site, nor_variants

_MAX_REMOVED_PER_CONNECTOR = 1


def _decompose(adj):
    """(junctions, edges): edges are (kind, u, v, atoms) with kind 'link' between junctions or 'term' toward a free end."""
    junction = {a for a, n in adj.items() if len(n) >= 3}
    if not junction:
        return None
    edges, seen = [], set()
    for u in sorted(junction, key=loc_key):
        for n in sorted(adj[u], key=loc_key):
            path, previous, current = [], u, n
            while current not in junction and len(adj[current]) == 2:
                path.append(current)
                previous, current = current, next(x for x in adj[current] if x != previous)
            if current in junction:
                key = frozenset(path) if path else frozenset((u, current))
                if key in seen:
                    continue
                seen.add(key)
                edges.append(("link", u, current, tuple(path)))
            else:
                edges.append(("term", u, None, tuple(path + [current])))
    return junction, edges


def _chains(view, start, length, used, end=None):
    """Simple paths of `length` unused molecule atoms leaving `start`; with `end`, the last atom must touch it."""
    if length == 0:
        if end is None or end in view.adj[start]:
            yield ()
        return
    stack = [(start, ())]
    while stack:
        atom, path = stack.pop()
        for n in view.adj[atom]:
            if n in used or n in path or n == start:
                continue
            nxt = path + (n,)
            if len(nxt) == length:
                if end is None or end in view.adj[n]:
                    yield nxt
            else:
                stack.append((n, nxt))


def _reaches(view, start, length, used):
    return length == 0 or next(_chains(view, start, length, used), None) is not None


def _reach(view, bonds):
    """{atom: {k: atoms joined to it by a simple path of k bonds}} for k up to `bonds`."""
    table = getattr(view, "reach", None)
    if table is not None and view.reach_bonds >= bonds:
        return table
    table = {}
    for start in view.adj:
        ends = {}
        stack = [(start, (start,))]
        while stack:
            atom, path = stack.pop()
            for n in view.adj[atom]:
                if n in path:
                    continue
                ends.setdefault(len(path), set()).add(n)
                if len(path) < bonds:
                    stack.append((n, path + (n,)))
        table[start] = ends
    view.reach, view.reach_bonds = table, bonds
    return table


def _least_change(table, adj, a, b, length, limit, floor=0):
    """Smallest |change| of a connector of `length` atoms that a simple path from a to b can realise, None if none.

    A connector that leaves a junction and returns to it (a == b) is a cycle through a of at least three atoms."""
    best = None
    for change in range(-_MAX_REMOVED_PER_CONNECTOR, limit + 1):
        atoms = length + change
        if atoms < floor or (best is not None and abs(change) >= best):
            continue
        if a == b:
            reachable = atoms >= 2 and any(end in adj[a] for end in table[a].get(atoms, ()))
        else:
            reachable = b in table[a].get(atoms + 1, ())
        if reachable:
            best = abs(change)
    return best


def _domains(junction, around, adj, view, table, budget, cuttable):
    """Arc-consistent candidate molecule atoms of every junction; None when one is empty."""
    lost = {
        a: max(
            1 if budget else 0,
            min(budget, sum(1 for n in adj[a] if len(adj[n]) == 1)) + (1 if any(i in cuttable for i, _, _, _ in around[a]) else 0),
        )
        for a in junction
    }
    needed = {a: max(1, len(adj[a]) - lost[a]) for a in junction}
    spans = {}

    def span(atom, length):
        if (atom, length) not in spans:
            ends = set()
            for change in range(-_MAX_REMOVED_PER_CONNECTOR, budget + 1):
                if length + change >= 0:
                    ends |= table[atom].get(length + change + 1, set())
            spans[(atom, length)] = ends
        return spans[(atom, length)]

    domains = {u: {m for m in view.adj if len(view.adj[m]) >= needed[u]} for u in junction}
    changed = True
    while changed:
        changed = False
        for u in junction:
            keep = {
                m for m in domains[u]
                if all(
                    span(m, n) & domains[b if a == u else a]
                    for i, a, b, n in around[u]
                    if (b if a == u else a) != u and i not in cuttable
                )
            }
            if not keep:
                return None
            if len(keep) != len(domains[u]):
                domains[u] = keep
                changed = True
    return domains


_MAX_SEED_MATCHES = 20000
WORK_LIMIT = 150000
WORK_PER_PARENT = 25000


class Exhausted(Exception):
    """The search effort allowed for one molecule is spent."""

_MIN_ATOMS_PER_OPERATION = 5
_MAX_COST = 4


def _on_cycle(adj, edge):
    """True when the connector lies on a ring (its end junctions stay joined without it)."""
    _, u, v, atoms = edge
    removed = set(atoms)
    if u == v:
        return True
    seen, stack = {u}, [u]
    while stack:
        a = stack.pop()
        for n in adj[a]:
            if n in removed or n in seen or (a == u and n == v and not atoms):
                continue
            if n == v:
                return True
            seen.add(n)
            stack.append(n)
    return False


def _order(edges):
    """Connector indices in breadth-first order over the junctions, so every prefix is a connected piece of the parent."""
    start = edges[0][1] if edges else None
    order, reached, queue = [], {start}, [start]
    remaining = list(range(len(edges)))
    while remaining:
        if not queue:
            queue.append(next(a for i in remaining for a in (edges[i][1], edges[i][2]) if a is not None))
            reached.add(queue[0])
        junction = queue.pop(0)
        for i in list(remaining):
            kind, u, v, atoms = edges[i]
            if junction in (u, v):
                order.append(i)
                remaining.remove(i)
                other = v if junction == u else u
                if other is not None and other not in reached:
                    reached.add(other)
                    queue.append(other)
    return order


def _runs(edges, order, count):
    """`count` consecutive runs of `order` that hold about the same number of atoms (some may be empty).

    Pigeonhole: k operations each damage one connector, so at least one of k + 1 runs stays intact."""
    total = sum(len(edges[i][3]) + 1 for i in order)
    groups, size = [[]], 0
    for i in order:
        if len(groups) < count and size >= total * len(groups) / count:
            groups.append([])
        groups[-1].append(i)
        size += len(edges[i][3]) + 1
    return groups + [[] for _ in range(count - len(groups))]


def _partition(edges, count):
    return _runs(edges, _order(edges), count)


def _group_query(edges, group):
    chains = [([u] + list(path) + ([v] if kind == "link" else [])) for kind, u, v, path in (edges[i] for i in group)]
    atoms = sorted({a for chain in chains for a in chain}, key=loc_key)
    index = {a: k for k, a in enumerate(atoms)}
    query = Chem.RWMol()
    for _ in atoms:
        query.AddAtom(Chem.AtomFromSmarts("*"))
    for bond in {frozenset((a, b)) for chain in chains for a, b in zip(chain, chain[1:])}:
        a, b = tuple(bond)
        query.AddBond(index[a], index[b], Chem.BondType.SINGLE)
    return atoms, index, query.GetMol()


_IMPOSSIBLE = 99


def _shape_bound(adj, view):
    """Operations the skeleton shape forces: none adds a bond, a nor removes one, the single seco removes two and a ring."""
    molecule = view.__dict__.get("degrees")
    if molecule is None:
        molecule = view.degrees = sorted((len(n) for n in view.adj.values()), reverse=True)
    parent = sorted((len(n) for n in adj.values() if len(n) >= 3), reverse=True)
    short = sum(max(0, d - (molecule[i] if i < len(molecule) else 0)) for i, d in enumerate(parent))
    rings = sum(len(n) for n in adj.values()) // 2 - len(adj) + 1
    if rings - 1 > view.cyclomatic:
        return _IMPOSSIBLE
    return max(0, short - (rings > 0))


def lower_bound(parent, view, limit):
    """Fewest operations that can still leave one of the pigeonhole groups intact (limit + 1 when none can)."""
    cache = view.__dict__.setdefault("lower", {})
    if parent.name not in cache:
        plan = _decompose(parent.adj)
        cache[parent.name] = 0 if plan is None else limit + 1
        if plan is not None:
            edges = plan[1]
            for k in range(max(1, _shape_bound(parent.adj, view)), limit + 1):
                if len(parent.adj) < _MIN_ATOMS_PER_OPERATION * (k + 1):
                    break
                groups = _partition(edges, k + 1)
                if any(not g for g in groups) or any(view.flat.HasSubstructMatch(_group_query(edges, g)[2]) for g in groups):
                    cache[parent.name] = k
                    break
    return cache[parent.name]


_MAX_CANDIDATE_TESTS = 400
_SEEDS_BEFORE_NARROWING = 2000


def _unmatched(edges, view, group, memo):
    key = frozenset(group)
    if key not in memo:
        memo[key] = not view.flat.HasSubstructMatch(_group_query(edges, sorted(key))[2])
    return memo[key]


def _components(edges, connectors):
    """Connected pieces of a set of connectors: two connectors are joined by a junction they share.

    The matcher is only given connected pieces; an unconnected query makes it try every combination of placements."""
    pieces, remaining = [], list(connectors)
    while remaining:
        piece = [remaining.pop(0)]
        junctions = {j for j in edges[piece[0]][1:3] if j is not None}
        grew = True
        while grew:
            grew = False
            for i in list(remaining):
                ends = {j for j in edges[i][1:3] if j is not None}
                if ends & junctions:
                    piece.append(i)
                    remaining.remove(i)
                    junctions |= ends
                    grew = True
        pieces.append(piece)
    return pieces


def _shrink(edges, view, group, memo):
    """An irreducible connected part of the connected `group` that the molecule still cannot contain unchanged."""
    group = list(group)
    shrunk = True
    while shrunk:
        shrunk = False
        for i in list(group):
            rest = [j for j in group if j != i]
            if rest and len(_components(edges, rest)) == 1 and _unmatched(edges, view, rest, memo):
                group, shrunk = rest, True
                break
    return frozenset(group)


def _obstructions(edges, view, memo):
    """Irreducible sets of connectors that the molecule cannot contain unchanged: an operation must change one of each."""
    touching, found = {}, set()
    for i, (_, u, v, _) in enumerate(edges):
        for junction in (u, v):
            if junction is not None:
                touching.setdefault(junction, []).append(i)
    for start in range(len(edges)):
        order, seen = [start], {start}
        for i in order:
            for junction in (edges[i][1], edges[i][2]):
                for nxt in touching.get(junction, ()):
                    if nxt not in seen:
                        seen.add(nxt)
                        order.append(nxt)
        if not _unmatched(edges, view, order, memo):
            continue
        low, high = 1, len(order)
        while low < high:
            middle = (low + high) // 2
            if _unmatched(edges, view, order[:middle], memo):
                high = middle
            else:
                low = middle + 1
        found.add(_shrink(edges, view, order[:low], memo))
    return found


def _hitting_sets(obstructions, limit):
    """Every set of at most `limit` connectors that meets all obstructions and is built from them, smallest first."""
    ordered = sorted(obstructions, key=lambda o: (len(o), sorted(o)))
    found = set()

    def grow(chosen):
        unhit = next((o for o in ordered if not o & chosen), None)
        if unhit is None:
            found.add(chosen)
        elif len(chosen) < limit:
            for i in sorted(unhit):
                grow(chosen | {i})

    grow(frozenset())
    return sorted(found, key=lambda s: (len(s), sorted(s)))


def _minimal_changes(edges, view, limit, key):
    """The smallest sets of connectors an operation may change so that the rest of the parent still fits unchanged, or None.

    A set is valid when every connected piece of the remainder fits, and every superset of a valid set is valid, so the
    minimal valid sets are the minimal hitting sets of all obstructions. Obstructions are learned lazily: a candidate
    with a piece that does not fit yields a new obstruction that rules out every candidate sharing that piece. Pieces
    are tested one by one, so the test is necessary but does not require them to avoid each other."""
    store = view.__dict__.setdefault("obstruction_store", {})
    if key is None or key not in store:
        memo = {}
        record = (memo, _obstructions(edges, view, memo))
        if key is not None:
            store[key] = record
    else:
        record = store[key]
    memo, family = record
    everything, tests = set(range(len(edges))), 0
    while True:
        valid = []
        for chosen in _hitting_sets(family, limit):
            if any(known <= chosen for known in valid):
                continue
            tests += 1
            if tests > _MAX_CANDIDATE_TESTS:
                return None
            failed = next((p for p in _components(edges, sorted(everything - chosen)) if _unmatched(edges, view, p, memo)), None)
            if failed is None:
                valid.append(chosen)
            else:
                family.add(_shrink(edges, view, failed, memo))
                break
        else:
            return valid




def _placements(edges, view, piece, cache):
    """(index, matches, query) of one connected piece; the matches are capped, which is enough to compare pieces."""
    key = tuple(sorted(piece))
    if key not in cache:
        _, index, query = _group_query(edges, key)
        cache[key] = (index, view.flat.GetSubstructMatches(query, uniquify=False, maxMatches=_SEEDS_BEFORE_NARROWING), query)
    return cache[key]


def _piece_seeds(edges, junction, view, piece, domains, cache, frozen):
    """Placements of the connectors of one connected piece left unchanged; the others stay open unless `frozen`."""
    index, matches, query = _placements(edges, view, piece, cache)
    if len(matches) >= _SEEDS_BEFORE_NARROWING:
        matches = view.flat.GetSubstructMatches(query, uniquify=False, maxMatches=_MAX_SEED_MATCHES)
    placed_at = [(u, index[u]) for u in junction if u in index]
    limits = [(position, domains[u]) for u, position in placed_at] if domains else []
    unchanged = {i: 0 for i in piece}
    opened = [i for i in range(len(edges)) if i not in unchanged]
    seeds, seen = [], set()
    for match in matches:
        if any(match[position] not in allowed for position, allowed in limits):
            continue
        mark = (tuple(match[position] for _, position in placed_at), frozenset(match))
        if mark not in seen:
            seen.add(mark)
            seeds.append(({u: match[position] for u, position in placed_at}, mark[1], unchanged, opened, frozen))
    return seeds


def _pigeonhole_seeds(edges, junction, view, budget):
    """Placements of one intact group out of budget + 1 equal slices of the parent, for when no obstruction guides the search."""
    everything = ({}, frozenset(), {}, list(range(len(edges))), None)
    groups = _partition(edges, budget + 1)
    if any(not group for group in groups):
        return [everything]
    seeds = []
    for group in groups:
        atoms, index, query = _group_query(edges, group)
        for match in view.flat.GetSubstructMatches(query, uniquify=False, maxMatches=_MAX_SEED_MATCHES):
            image = {a: match[index[a]] for a in atoms}
            seeds.append(
                (
                    {u: image[u] for u in junction if u in image},
                    frozenset(image.values()),
                    {i: 0 for i in group},
                    [i for i in range(len(edges)) if i not in group],
                    None,
                )
            )
    return seeds


def _seeds(edges, junction, view, budget, key=None, domains=None):
    """Starting states of the search, from connected pieces of the parent that cannot have changed.

    A solution changes a set of connectors that holds a minimal one M, and at most r = budget - |M| more. Pieces of the
    rest are chosen, fewest placements first, until no r connectors can touch all of them: one piece is then intact,
    and it is matched exactly. With r = 0 the solution changes exactly M, so every other connector is frozen."""
    minimal = _minimal_changes(edges, view, budget, key)
    if minimal is None:
        return _pigeonhole_seeds(edges, junction, view, budget)
    everything = frozenset(range(len(edges)))
    seeds, chosen, cache = [], set(), {}
    for changed in minimal:
        remainder = everything - changed
        frozen = remainder if len(changed) == budget else None
        pieces = []
        while True:
            touching = next(iter(_hitting_sets(pieces, budget - len(changed))), None)
            if touching is None:
                break
            candidates = _components(edges, sorted(remainder - touching))
            if not candidates:
                return [({}, frozenset(), {}, sorted(everything), None)]
            pieces.append(frozenset(min(candidates, key=lambda piece: len(_placements(edges, view, piece, cache)[1]))))
        for piece in pieces:
            if (piece, frozen) not in chosen:
                chosen.add((piece, frozen))
                seeds.extend(_piece_seeds(edges, junction, view, piece, domains, cache, frozen))
    return seeds



def _search(adj, view, budget, cuts=True, key=None):
    """(edges, deltas): the atom-count change of every connector, for each way of costing exactly `budget`."""
    plan = _decompose(adj)
    if plan is None:
        return None
    junction, edges = plan
    found = set()
    links = [(i, u, v, len(atoms)) for i, (kind, u, v, atoms) in enumerate(edges) if kind == "link"]
    table = _reach(view, max([len(atoms) for _, _, _, atoms in edges] + [1]) + budget + 2)
    around = {a: [(i, u, v, n) for i, u, v, n in links if a in (u, v)] for a in junction}
    ring_edge = {i: _on_cycle(adj, edges[i]) for i in range(len(edges)) if edges[i][0] == "link"}
    cuttable = {i for i, ring in ring_edge.items() if ring} if cuts and budget >= 1 else set()
    domains = _domains(junction, around, adj, view, table, budget, cuttable)
    if domains is None:
        return edges, found
    low = {i: 2 if edges[i][1] == edges[i][2] else 0 for i in range(len(edges))}
    visited = set()

    fixed = [None]

    def edge_options(i, phi, taken, room, cut_done):
        kind, u, v, atoms = edges[i]
        options = []
        frozen = fixed[0] is not None and i in fixed[0]
        for change in range(0 if frozen else -_MAX_REMOVED_PER_CONNECTOR, 1 if frozen else room + 1):
            new = len(atoms) + change
            if new < low[i]:
                continue
            if kind == "term":
                if _reaches(view, phi[u], new, taken):
                    options.append((change, abs(change), ()))
            else:
                options.extend((change, abs(change), path) for path in _chains(view, phi[u], new, taken, phi[v]))
        if cuts and kind == "link" and ring_edge[i] and room >= 1 and not cut_done and not frozen:
            for j in range(len(atoms) + 1):
                if _reaches(view, phi[u], j, taken) and _reaches(view, phi[v], len(atoms) - j, taken):
                    options.append((("cut", j), 1, ()))
        return options

    def junction_options(u, phi, taken, room):
        options = []
        for m in domains[u] - taken:
            extra = 0
            for i, a, b, n in around[u]:
                w = b if a == u else a
                if w != u and w in phi:
                    least = _least_change(table, view.adj, m, phi[w], n, room - extra)
                    if i in cuttable and (least is None or least > 1):
                        least = 1
                    if least is None:
                        extra = room + 1
                        break
                    extra += least
            if extra <= room:
                options.append(m)
        return options

    spent = [0]

    def solve(phi, taken, cost, deltas, open_edges):
        view.work -= 1
        spent[0] += 1
        if view.work < 0 or spent[0] > WORK_PER_PARENT:
            raise Exhausted()
        key = (tuple(sorted(deltas.items())), taken, frozenset(phi.items()), id(fixed[0]))
        if key in visited:
            return
        visited.add(key)
        room = budget - cost
        ready = [i for i in open_edges if edges[i][1] in phi and (edges[i][2] is None or edges[i][2] in phi)]
        if ready:
            cut_done = any(isinstance(d, tuple) for d in deltas.values())
            lower = {}
            for i in ready:
                if edges[i][0] != "link":
                    continue
                least = _least_change(table, view.adj, phi[edges[i][1]], phi[edges[i][2]], len(edges[i][3]), room, low[i])
                if fixed[0] is not None and i in fixed[0]:
                    least = 0 if least == 0 else None
                elif cuts and ring_edge[i] and not cut_done and (least is None or least > 1):
                    least = 1
                lower[i] = least
            if any(least is None for least in lower.values()) or sum(lower.values()) > room:
                return
            chosen = None
            for i in ready:
                options = edge_options(i, phi, taken, room - sum(lower.values()) + (lower.get(i) or 0), cut_done)
                if not options:
                    return
                if chosen is None or len(options) < len(chosen[1]):
                    chosen = (i, options)
            i, options = chosen
            rest = [j for j in open_edges if j != i]
            for delta, price, path in options:
                solve(phi, taken | frozenset(path), cost + price, {**deltas, i: delta}, rest)
            return
        pending = [u for u in junction if u not in phi]
        if not pending:
            if cost == budget and not open_edges:
                found.add(tuple(sorted(deltas.items())))
            return
        chosen = None
        for u in pending:
            options = junction_options(u, phi, taken, room)
            if not options:
                return
            if chosen is None or len(options) < len(chosen[1]):
                chosen = (u, options)
        u, options = chosen
        for m in options:
            solve({**phi, u: m}, taken | {m}, cost, deltas, open_edges)

    for phi, taken, deltas, rest, frozen in _seeds(edges, junction, view, budget, key, domains):
        if all(m in domains[u] for u, m in phi.items()):
            fixed[0] = frozen
            solve(phi, taken, 0, deltas, rest)
    return edges, found


def _sites_for(skel, atoms, ends, sites, nors):
    touched = set(atoms) | set(ends)
    homo = [site for site in sites if (({site[1]} if site[0] == "terminal" else set(site[1])) <= touched) and _site_alive(skel, site)]
    removable = [a for a in nors if a in atoms and a in skel.adj]
    return homo, removable


def _apply(skel, move):
    built = skel
    kind, target, count = move
    for _ in range(count):
        if kind == "homo":
            if not _site_alive(built, target):
                return None
            built = apply_homo(built, target)
        else:
            if target not in built.adj or len(built.adj[target]) > 2:
                return None
            built = apply_nor(built, target)
    return built


def _skeletons(skel, edges, results, sites, nors, produced):
    for deltas in results:
        moves, start = [], skel
        for i, change in deltas:
            if change == 0:
                continue
            _, u, v, atoms = edges[i]
            if isinstance(change, tuple):
                chain = [u] + list(atoms) + [v]
                start = apply_seco(start, *sorted((chain[change[1]], chain[change[1] + 1]), key=loc_key))
                continue
            homo, removable = _sites_for(skel, atoms, [u] if v is None else [u, v], sites, nors)
            moves.append([("homo", s, change) for s in homo] if change > 0 else [("nor", a, 1) for a in removable])
        if any(not options for options in moves):
            continue
        stack = [start]
        for options in moves:
            stack = [built for current in stack for option in options if (built := _apply(current, option)) is not None]
        for result in stack:
            produced.setdefault(len(result.ops), {})[tuple(sorted(str(o[:3]) for o in result.ops))] = result


def read_operations(parent, view, cost, terminal_only=False):
    """Skeletons with exactly `cost` nor/homo/seco operations that may explain the molecule; None without a junction."""
    if len(parent.adj) - cost > len(view.adj) or len(parent.adj) < _MIN_ATOMS_PER_OPERATION * (cost + 1):
        return []
    if cost < lower_bound(parent, view, _MAX_COST) or view.work < 0:
        return []
    sites = [site for site in homo_sites(parent) if not is_branch_site(parent, site)]
    nors = nor_variants(parent)
    if terminal_only:
        sites = [s for s in sites if s[0] == "terminal"]
        nors = [a for a in nors if len(parent.adj[a]) == 1]
    skel = Skel.of(parent)
    try:
        outcome = _search(skel.adj, view, cost, not terminal_only, parent.name)
    except Exhausted:
        return []
    if outcome is None:
        return None
    edges, results = outcome
    produced = {}
    _skeletons(skel, edges, results, sites, nors, produced)
    return list(produced.get(cost, {}).values())
