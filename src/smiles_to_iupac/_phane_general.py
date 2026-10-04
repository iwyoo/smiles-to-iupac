"""Phane parent hydrides of benzene amplificants joined by single-atom or chain bridges, linear or cyclic, with
substituents on the rings (composite locants 14, P-26.4.3) and on the bridge atoms: 1,4(1,4)-dibenzenacyclohexaphane,
12-bromo-1,4(1,4)-dibenzenacyclohexaphane, 2,4,6-trioxa-1,7(1),3,5(1,4)-tetrabenzenaheptaphane (P-26, P-52.2.5.1)."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._fusion_numbering_general import _HETERO_RANK
from ._numerals import multiplying_prefix, numerical_term
from ._substituents import format_substituent_prefixes, name_branch

_PREFIX = {8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza", 14: "sila", 15: "phospha", 33: "arsa", 32: "germa", 5: "bora"}


class PhaneLoc(int):
    """Locant of a phane skeleton atom or, with `local`, of an amplificant atom (primary 1, local 4 -> '14')."""

    def __new__(cls, primary, local=0):
        obj = super().__new__(cls, primary * 1000 + local)
        obj.primary, obj.local = primary, local
        return obj

    def __str__(self):
        return f"{self.primary}{self.local}" if self.local else str(self.primary)

    __repr__ = __str__
    __format__ = lambda self, spec: str(self)


def _benzene_rings(mol):
    return [
        list(ring)
        for ring in mol.GetRingInfo().AtomRings()
        if len(ring) == 6 and all(mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring)
    ]


def _ring_order(graph, ring):
    cycle = [ring[0]]
    previous = None
    while len(cycle) < 6:
        following = next(n for n in graph[cycle[-1]] if n in ring and n != previous and n not in cycle)
        previous = cycle[-1]
        cycle.append(following)
    return cycle


def _component(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n not in blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _path(graph, source, target, allowed):
    parent = {source: None}
    queue = [source]
    while queue:
        a = queue.pop(0)
        if a == target:
            break
        for n in graph[a]:
            if n in allowed and n not in parent:
                parent[n] = a
                queue.append(n)
    if target not in parent:
        return None
    path = [target]
    while path[-1] != source:
        path.append(parent[path[-1]])
    return path[::-1]


def find_phane(mol):
    """(rings, bridges, substituent roots, cyclic) of a supported phane, else None. A bridge is
    (ring_a, atom_a, [skeleton atoms], atom_b, ring_b)."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    rings = _benzene_rings(mol)
    if rings is None or len(rings) < 2:
        return None
    ring_atoms = {a for ring in rings for a in ring}
    if len(ring_atoms) != 6 * len(rings):
        return None
    ring_of = {a: i for i, ring in enumerate(rings) for a in ring}
    graph = adjacency(mol)
    seen, bridges, branches = set(), [], []
    for a in sorted(ring_atoms):
        for n in sorted(graph[a]):
            if n in ring_atoms and ring_of[n] != ring_of[a] and (n, a) not in seen:
                seen.add((a, n))
                bridges.append((ring_of[a], a, [], n, ring_of[n]))
    chain = [a for a in range(mol.GetNumAtoms()) if a not in ring_atoms]
    claimed = set()
    for start in chain:
        if start in claimed:
            continue
        component = _component(graph, start, ring_atoms)
        claimed |= component
        contacts = [(r, c) for c in component for r in graph[c] if r in ring_atoms]
        if len(contacts) == 1:
            branches.append((contacts[0][0], contacts[0][1], component))
            continue
        if len(contacts) != 2 or ring_of[contacts[0][0]] == ring_of[contacts[1][0]]:
            return None
        (ra, ca), (rb, cb) = contacts
        path = _path(graph, ca, cb, component)
        if path is None or any(mol.GetBondBetweenAtoms(x, y).GetBondTypeAsDouble() != 1.0 for x, y in zip(path, path[1:])):
            return None
        bridges.append((ring_of[ra], ra, path, rb, ring_of[rb]))
        skeleton = set(path)
        for atom in path:
            for n in graph[atom]:
                if n in component and n not in skeleton:
                    sub = _component(graph, n, skeleton | ring_atoms)
                    if any(x in sub for x in sum(([c for c in graph[s] if c in ring_atoms] for s in sub), [])):
                        return None
                    branches.append((atom, n, sub))
    degree = {i: 0 for i in range(len(rings))}
    for ring_a, _, _, _, ring_b in bridges:
        degree[ring_a] += 1
        degree[ring_b] += 1
    if any(d not in (1, 2) for d in degree.values()):
        return None
    ends = [i for i, d in degree.items() if d == 1]
    if len(ends) not in (0, 2) or len(bridges) != len(rings) - (1 if ends else 0):
        return None
    cyclic = not ends
    if cyclic:
        for index, ring in enumerate(rings):
            attached = [a for ra, a, _, _, rb in bridges if ra == index] + [b for ra, _, _, b, rb in bridges if rb == index]
            cycle = _ring_order(graph, ring)
            gap = abs(cycle.index(attached[0]) - cycle.index(attached[1]))
            if min(gap, 6 - gap) < 2:
                return None
    return rings, bridges, branches, cyclic


def _walks(rings, bridges, cyclic):
    """Every skeleton traversal as a node list: ('ring', i, entry atom, exit atom) and ('atom', atom)."""
    adjacency_of = {i: [] for i in range(len(rings))}
    for k, (ra, _, _, _, rb) in enumerate(bridges):
        adjacency_of[ra].append((k, False))
        adjacency_of[rb].append((k, True))
    starts = list(adjacency_of) if cyclic else [i for i in adjacency_of if len(adjacency_of[i]) == 1]
    results = []
    for start in starts:
        for first in range(len(adjacency_of[start])):
            nodes, ring, used, entry = [], start, set(), None
            while True:
                options = [o for o in adjacency_of[ring] if o[0] not in used]
                if ring == start and not nodes:
                    options = [adjacency_of[ring][first]]
                if not options:
                    nodes.append(("ring", ring, entry, None))
                    break
                k, reverse = options[0]
                ra, aa, path, ab, rb = bridges[k]
                exit_atom, chain, entry_next, next_ring = (ab, path[::-1], aa, ra) if reverse else (aa, path, ab, rb)
                nodes.append(("ring", ring, entry, exit_atom))
                nodes.extend(("atom", a) for a in chain)
                used.add(k)
                entry = entry_next
                if cyclic and next_ring == start:
                    nodes[0] = ("ring", start, entry, nodes[0][3])
                    break
                ring = next_ring
            results.append(nodes)
    return results


def _ring_numbering(graph, ring, attachments, lower_attachment, substituent_atoms, free_atom=None):
    """Best local numbering of a benzene amplificant: {atom: local locant} and the attachment locants cited."""
    cycle = _ring_order(graph, ring)
    best = None
    for start in range(6):
        for step in (1, -1):
            order = [cycle[(start + step * k) % 6] for k in range(6)]
            local = {a: i + 1 for i, a in enumerate(order)}
            att = sorted(local[a] for a in attachments)
            subs = sorted(local[a] for a in substituent_atoms)
            adjacency_rule = 0 if lower_attachment is None or len(attachments) < 2 or local[lower_attachment] < min(local[a] for a in attachments if a != lower_attachment) else 1
            free = [local[free_atom]] if free_atom in local else []
            key = (att, free, subs, adjacency_rule)
            if best is None or key < best[0]:
                best = (key, local)
    return best[1]


def _evaluate(mol, graph, rings, bridges, branches, nodes, cyclic, halogens, free_atom=None):
    nodes_flat = nodes
    size = len(nodes_flat)
    position = {}
    ring_position = {}
    for index, node in enumerate(nodes_flat):
        if node[0] == "atom":
            position[node[1]] = index + 1
        else:
            ring_position[node[1]] = index + 1
    ring_sub_atoms = {i: [] for i in range(len(rings))}
    chain_subs = []
    free_loc = None
    ring_of = {a: i for i, ring in enumerate(rings) for a in ring}
    for owner, root, _ in branches:
        if owner in ring_of:
            ring_sub_atoms[ring_of[owner]].append((owner, root))
        else:
            chain_subs.append((position[owner], owner, root))
    composite = []
    locals_by_ring = {}
    attachment_text = {}
    for index, node in enumerate(nodes_flat):
        if node[0] != "ring":
            continue
        _, ring_index, entry, exit_atom = node
        attachments = [a for a in (entry, exit_atom) if a is not None]
        neighbours = []
        if entry is not None:
            neighbours.append((index if index > 0 else size, entry))
        if exit_atom is not None:
            neighbours.append((index + 2 if index + 2 <= size else 1, exit_atom))
        lower = min(neighbours)[1] if len(neighbours) == 2 else None
        local = _ring_numbering(graph, rings[ring_index], attachments, lower, [a for a, _ in ring_sub_atoms[ring_index]], free_atom)
        locals_by_ring[ring_index] = local
        attachment_text[index + 1] = tuple(sorted(local[a] for a in attachments))
        if free_atom in local:
            free_loc = PhaneLoc(index + 1, local[free_atom])
        for owner, root in ring_sub_atoms[ring_index]:
            composite.append((PhaneLoc(index + 1, local[owner]), owner, root))
    if free_atom in position:
        free_loc = PhaneLoc(position[free_atom])
    chain_locs = [(PhaneLoc(p), owner, root) for p, owner, root in chain_subs]
    all_subs = composite + chain_locs
    hetero = sorted(
        (_HETERO_RANK.get(mol.GetAtomWithIdx(node[1]).GetSymbol(), 99), index + 1)
        for index, node in enumerate(nodes_flat)
        if node[0] == "atom" and mol.GetAtomWithIdx(node[1]).GetAtomicNum() != 6
    )
    superatoms = sorted(ring_position.values())
    key = (
        superatoms,
        sorted(p for _, p in hetero),
        [p for _, p in sorted(hetero)],
        [attachment_text[p] for p in superatoms],
        [(free_loc.primary, int(free_loc))] if free_atom is not None else [],
        sorted(loc.primary for loc, _, _ in all_subs),
        sorted(int(loc) for loc, _, _ in all_subs),
    )
    return key, attachment_text, all_subs, (free_loc if free_atom is not None else None)


def name_phane_general(mol, free_atom=None):
    found = find_phane(mol)
    if found is None:
        raise UnsupportedStructure("this structure is not a supported phane")
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        raise UnsupportedStructure("stereo on a phane is not supported yet")
    rings, bridges, branches, cyclic = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    best = None
    for nodes in _walks(rings, bridges, cyclic):
        key, attachment_text, subs, free_loc = _evaluate(mol, graph, rings, bridges, branches, nodes, cyclic, halogens, free_atom)
        if best is None or key < best[0]:
            best = (key, nodes, attachment_text, subs, free_loc)
    _, nodes, attachment_text, subs, free_loc = best
    size = len(nodes)

    grouped = {}
    for loc, owner, root in subs:
        name, compound = name_branch(graph, root, owner, halogens, frozenset(), mol=mol, unsaturated=True)
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(loc)
    prefix = format_substituent_prefixes(grouped) if grouped else ""

    patterns = {}
    for position, att in attachment_text.items():
        patterns.setdefault(att, []).append(position)
    groups = ",".join(
        f"{','.join(map(str, sorted(locs)))}({','.join(map(str, att))})"
        for att, locs in sorted(patterns.items(), key=lambda kv: kv[0])
    )
    hetero_text = {}
    for index, node in enumerate(nodes):
        if node[0] == "atom":
            z = mol.GetAtomWithIdx(node[1]).GetAtomicNum()
            if z != 6:
                if z not in _PREFIX:
                    raise UnsupportedStructure("an unsupported replacement atom in a phane bridge")
                hetero_text.setdefault(z, []).append(index + 1)
    replacement = "-".join(
        f"{','.join(map(str, locs))}-{multiplying_prefix(len(locs)) if len(locs) > 1 else ''}{_PREFIX[z]}"
        for z, locs in sorted(hetero_text.items(), key=lambda kv: _HETERO_RANK[Chem.GetPeriodicTable().GetElementSymbol(kv[0])])
    )
    count = sum(len(v) for v in patterns.values())
    amplificant = f"{multiplying_prefix(count)}benzena"
    parent = (f"cyclo{numerical_term(size)}" if cyclic else numerical_term(size)) + "phane"
    core = (replacement + "-" if replacement else "") + f"{groups}-{amplificant}{parent}"
    if free_loc is not None:
        core = f"{core[:-1]}-{free_loc}-yl"
    return f"{prefix}-{core}" if prefix else core


def phane_substituent(mol, graph, root, coming_from):
    """(name, True) of a phane group entered at `root`, or None when the branch is not a supported phane."""
    arm = _component(graph, root, {coming_from})
    benzene = [r for r in _benzene_rings(mol) if set(r) <= arm]
    if len(benzene) < 2:
        return None
    editable = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - arm, reverse=True):
        editable.RemoveAtom(index)
    kept = sorted(arm)
    fragment = editable.GetMol()
    try:
        Chem.SanitizeMol(fragment)
    except Exception:
        return None
    free_atom = kept.index(root)
    found = find_phane(fragment)
    if found is None:
        return None
    rings, bridges, _, cyclic = found
    if not cyclic and (len(rings) < 4 or len(_walks(rings, bridges, cyclic)[0]) < 7):
        return None
    return name_phane_general(fragment, free_atom), True
