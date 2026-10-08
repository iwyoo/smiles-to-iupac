"""Phane parent hydrides (P-26): any ring or ring system as an amplificant, joined by atoms or chains into an
unbranched, monocyclic, von Baeyer or spiro simplified skeleton, with composite locants (P-26.4.3), skeletal
replacement in skeleton and amplificants (P-26.5), substituents, a principal group suffix and free valences."""

from contextvars import ContextVar
from itertools import combinations

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, elides_before, halogen_substituents, multiplied_word, sanitize_probe, specified_stereo_elements
from ._multiplicative_groups import SUFFIX_RANKS
from ._numerals import multiplying_prefix
from ._phane_amplificant import PhaneLoc, build_amplificant, needs_bis
from ._phane_skeleton import skeleton_numberings
from ._pin import mark
from ._ring_diyl_numbering import _PREFIX, _RANK
from ._substituents import BRANCH_STEREO, format_substituent_prefixes, name_branch

_LARGEST_AMPLIFICANT_RING_ALWAYS = 8
_MAX_LARGE_RINGS = 12
_CACHE = {}
_BUILDING = ContextVar("phane_building", default=False)


class Phane:
    """A decomposition of a molecule into amplificants and skeleton atoms: `nodes` are ('amp', index) or ('atom',
    atom); `attachments[i]` lists (atom of node i, other node) per skeleton bond of an amplificant."""

    def __init__(self, amplificants, nodes, attachments, branches, cyclic, skeleton, fused, multiple):
        self.amplificants = amplificants
        self.nodes = nodes
        self.attachments = attachments
        self.branches = branches
        self.cyclic = cyclic
        self.skeleton = skeleton
        self.fused = fused
        self.multiple = multiple
        self.size = len(nodes)


def is_pin_phane(phane):
    """P-52.2.5.1: a cyclophane has six nodes and a mancude amplificant, none of them ortho or peri fused to the
    skeleton ring; a linear phane has four amplificants and seven nodes."""
    if phane.cyclic:
        return not phane.fused and phane.size >= 6 and any(a.mancude for a in phane.amplificants)
    return len(phane.amplificants) >= 4 and phane.size >= 7


def is_fused_cyclophane(phane):
    """A cyclophane whose only departure from P-52.2.5.1 is ortho or peri fusion to the skeleton ring: the phane name is
    valid but a fusion or bridged fused name is preferred."""
    return phane.cyclic and phane.fused and phane.size >= 6 and any(a.mancude for a in phane.amplificants)


def _free_atoms(free_atom):
    if free_atom is None:
        return ()
    return tuple(free_atom) if isinstance(free_atom, (tuple, list)) else (free_atom,)


def _component(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n not in blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _prefix_only(mol, atoms):
    """True when a substituent group carries only characteristic groups that are cited as prefixes: a phane parent
    with a principal group (hydroxy, amino, oxo, carboxy, cyano, ...) takes it as a suffix (P-41, P-44.1.1)."""
    for a in atoms:
        atom = mol.GetAtomWithIdx(a)
        z = atom.GetAtomicNum()
        if z in (6, 9, 17, 35, 53):
            if z == 6 and any(
                b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16) for b in atom.GetBonds()
            ):
                return False
            continue
        if z in (8, 16) and atom.GetDegree() == 2 and atom.GetTotalNumHs() == 0:
            continue
        return False
    return True


_SUFFIX_WORDS = {
    "carboxylic_acid": "carboxylic acid",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
    "ketone": "one",
    "alcohol": "ol",
    "thiol": "thiol",
    "amine": "amine",
}


def _suffix_class(mol, owner, root, atoms):
    """Class of a branch that is exactly one principal characteristic group attached directly to the phane (P-41),
    else None."""
    atom = mol.GetAtomWithIdx(root)
    bond_order = mol.GetBondBetweenAtoms(owner, root).GetBondTypeAsDouble()
    z = atom.GetAtomicNum()
    if len(atoms) == 1 and atom.GetDegree() == 1:
        if bond_order == 2.0:
            return "ketone" if z == 8 and not mol.GetAtomWithIdx(owner).GetIsAromatic() else None
        if bond_order != 1.0:
            return None
        if z == 8 and atom.GetTotalNumHs() == 1:
            return "alcohol"
        if z == 16 and atom.GetTotalNumHs() == 1:
            return "thiol"
        if z == 7 and atom.GetTotalNumHs() == 2:
            return "amine"
        return None
    if z != 6 or bond_order != 1.0:
        return None
    others = [(n.GetAtomicNum(), mol.GetBondBetweenAtoms(root, n.GetIdx()).GetBondTypeAsDouble(), n.GetTotalNumHs(), n.GetDegree()) for n in atom.GetNeighbors() if n.GetIdx() != owner]
    if len(atoms) == 2 and others == [(7, 3.0, 0, 1)]:
        return "nitrile"
    if len(atoms) == 2 and others == [(8, 2.0, 0, 1)] and atom.GetTotalNumHs() == 1:
        return "aldehyde"
    if len(atoms) == 3 and sorted(others) == [(8, 1.0, 1, 1), (8, 2.0, 0, 1)]:
        return "carboxylic_acid"
    return None


def _principal_suffix(mol, branches):
    """(class, roots) of the senior suffix-capable class among the branches, or (None, frozenset())."""
    classes = {}
    for owner, root, atoms in branches:
        cls = _suffix_class(mol, owner, root, atoms)
        if cls is not None:
            classes.setdefault(cls, set()).add(root)
    if not classes:
        return None, frozenset()
    best = min(classes, key=lambda c: (SUFFIX_RANKS[c], c == "thiol"))
    return best, frozenset(classes[best])


def _merge_rings(subset, ring_atoms):
    """Rings of `subset` joined when they share an atom: [(atoms, ring indices)]."""
    groups = []
    for index in subset:
        atoms = set(ring_atoms[index])
        joined = [g for g in groups if g[0] & atoms]
        for g in joined:
            groups.remove(g)
            atoms |= g[0]
        groups.append((atoms, [index] + [i for g in joined for i in g[1]]))
    return groups


def _leaf_prune(adj, alive, amplificant_count, amplificants_too):
    """Remove nodes with at most one live neighbour: acyclic chain atoms, and amplificant nodes (pendant rings that
    become substituents) only when `amplificants_too`."""
    stack = list(alive)
    while stack:
        n = stack.pop()
        if n not in alive or sum(1 for x in adj[n] if x in alive) > 1 or (n < amplificant_count and not amplificants_too):
            continue
        alive.discard(n)
        stack.extend(x for x in adj[n] if x in alive)


def _attached_to_skeleton_ring(component_rings, atoms, attached, graph, joined):
    """True when two attachment atoms that close a skeleton ring (`joined(a, b)`) are ortho or peri positions: the
    skeleton ring would then be fused to the amplificant. Bridged ring systems have no peri positions."""
    ortho_fused = all(len(r & q) <= 2 for r, q in combinations(component_rings, 2))
    for (a, other_a), (b, other_b) in combinations(attached, 2):
        if not joined(other_a, other_b):
            continue
        if b in graph[a]:
            return True
        if ortho_fused and any(sum(1 for ring in component_rings if c in ring) >= 2 for c in (graph[a] & graph[b]) & atoms):
            return True
    return False


def _connected_without(count, edge_list, removed):
    parent = list(range(count))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j, _, _ in edge_list:
        if removed not in (i, j):
            parent[find(i)] = find(j)
    return lambda a, b: find(a) == find(b)


def _try_subset(mol, graph, subset, ring_atoms, ring_bonds, ignore, free_atom, free_order, known):
    groups = _merge_rings(subset, ring_atoms)
    comp_of = {a: i for i, (atoms, _) in enumerate(groups) for a in atoms}
    induced = [set() for _ in groups]
    for bond in mol.GetBonds():
        i, j = comp_of.get(bond.GetBeginAtomIdx()), comp_of.get(bond.GetEndAtomIdx())
        if i is not None and i == j:
            induced[i].add(bond.GetIdx())
    for i, (_, indices) in enumerate(groups):
        if induced[i] != set().union(*(ring_bonds[r] for r in indices)):
            return None
    count = len(groups)
    chain_atoms = [a for a in range(mol.GetNumAtoms()) if a not in comp_of and a not in ignore]
    node_of = dict(comp_of)
    node_of.update({a: count + k for k, a in enumerate(chain_atoms)})
    adj = {n: {} for n in range(count + len(chain_atoms))}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in ignore or b in ignore or node_of[a] == node_of[b]:
            continue
        u, v = node_of[a], node_of[b]
        if v in adj[u]:
            return None
        adj[u][v] = (a, b, bond)
        adj[v][u] = (b, a, bond)
    alive = set(adj)
    _leaf_prune(adj, alive, count, False)
    edges = sum(1 for n in alive for x in adj[n] if x in alive) // 2
    cyclic = edges - len(alive) + 1 >= 1
    if cyclic:
        _leaf_prune(adj, alive, count, True)
    elif sum(1 for n in alive if sum(1 for x in adj[n] if x in alive) == 1) != 2:
        return None
    core = sorted(alive)
    if sum(1 for n in core if n < count) < 2:
        return None
    index = {n: i for i, n in enumerate(core)}
    edge_list, multiple = [], []
    for n in core:
        for x, (a, b, bond) in adj[n].items():
            if x in alive and n < x:
                if bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                    return None
                edge_list.append((index[n], index[x], a, b))
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    multiple.append((index[n], a, index[x], b, int(bond.GetBondTypeAsDouble()), bond.GetIdx()))
    core_atoms = set()
    for n in core:
        core_atoms |= groups[n][0] if n < count else {chain_atoms[n - count]}
    if any(a not in core_atoms for a in _free_atoms(free_atom)):
        return None
    branches = []
    for n in core:
        owner_atoms = groups[n][0] if n < count else {chain_atoms[n - count]}
        for a in owner_atoms:
            for r in graph[a]:
                if r in ignore or r in core_atoms:
                    continue
                branches.append((a, r, _component(graph, r, core_atoms | ignore)))
    if any(
        not _prefix_only(mol, atoms) and (free_atom is not None or _suffix_class(mol, owner, root, atoms) is None)
        for owner, root, atoms in branches
    ):
        return None
    attachments = {index[n]: [] for n in core if n < count}
    for i, j, a, b in edge_list:
        if i in attachments:
            attachments[i].append((a, j))
        if j in attachments:
            attachments[j].append((b, i))
    fused = cyclic and any(
        n < count
        and _attached_to_skeleton_ring(
            [ring_atoms[r] for r in groups[n][1]],
            groups[n][0],
            attachments[index[n]],
            graph,
            _connected_without(len(core), edge_list, index[n]),
        )
        for n in core
    )
    amplificants, nodes = [], []
    for n in core:
        if n >= count:
            nodes.append(("atom", chain_atoms[n - count]))
            continue
        key = frozenset(groups[n][0])
        if key not in known:
            known[key] = build_amplificant(mol, key, free_atom, free_order)
        if known[key] is None:
            return None
        nodes.append(("amp", len(amplificants)))
        amplificants.append(known[key])
    for i, _, j, _, _, _ in multiple:
        if any(nodes[k][0] == "amp" and not amplificants[nodes[k][1]].saturated_parent for k in (i, j)):
            return None
    skeleton = skeleton_numberings(len(core), [(i, j) for i, j, _, _ in edge_list])
    if not skeleton:
        return None
    return Phane(amplificants, nodes, attachments, branches, cyclic, skeleton, fused, multiple)


def _signature(mol, free_atom, free_order, ignore):
    atoms = tuple((a.GetAtomicNum(), a.GetIsAromatic(), a.GetTotalNumHs()) for a in mol.GetAtoms())
    bonds = tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType())) for b in mol.GetBonds())
    return atoms, bonds, _free_atoms(free_atom), free_order, tuple(sorted(ignore))


def find_phane(mol, free_atom=None, free_order=1, ignore=frozenset()):
    """The Phane decomposition of `mol` with the most amplificant rings, else None."""
    if _BUILDING.get():
        return None
    key = _signature(mol, free_atom, free_order, ignore)
    if key in _CACHE:
        return _CACHE[key]
    found = _search(mol, free_atom, free_order, ignore)
    if len(_CACHE) > 256:
        _CACHE.clear()
    _CACHE[key] = found
    return found


def _search(mol, free_atom, free_order, ignore):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    graph = {a: {n for n in ns if n not in ignore} for a, ns in graph.items() if a not in ignore}
    graph.update({a: set() for a in ignore})
    atom_count = sum(1 for a in range(mol.GetNumAtoms()) if a not in ignore)
    bond_count = sum(len(ns) for ns in graph.values()) // 2
    if bond_count - atom_count + 1 < 3:
        return None
    info = mol.GetRingInfo()
    ring_atoms = [frozenset(r) for r in info.AtomRings()]
    ring_bonds = [frozenset(r) for r in info.BondRings()]
    small = tuple(i for i, r in enumerate(ring_atoms) if len(r) <= _LARGEST_AMPLIFICANT_RING_ALWAYS)
    large = [i for i, r in enumerate(ring_atoms) if len(r) > _LARGEST_AMPLIFICANT_RING_ALWAYS]
    if not small or len(large) > _MAX_LARGE_RINGS:
        return None
    known = {}
    fused = None
    token = _BUILDING.set(True)
    try:
        for size in range(len(large), -1, -1):
            parents = []
            for chosen in combinations(large, size):
                found = _try_subset(mol, graph, small + chosen, ring_atoms, ring_bonds, ignore, free_atom, free_order, known)
                if found is not None and not found.fused:
                    parents.append(found)
                fused = fused or found
            if len(parents) > 1:
                return min(parents, key=lambda p: phane_seniority_key(mol, p) or ())
            if parents:
                return parents[0]
    finally:
        _BUILDING.reset(token)
    return fused


_CODE_RANK = {"Z": 0, "E": 1, "R": 0, "S": 1, "r": 0, "s": 1}


def _hetero_key(entries):
    """Heteroatom locants without regard to kind, then by seniority of the element (P-26.5.4.1, P-26.5.4.2)."""
    primary = sorted(loc.primary for _, loc in entries)
    complete = sorted(int(loc) for _, loc in entries)
    by_element = []
    for rank in sorted({r for r, _ in entries}):
        locs = [loc for r, loc in entries if r == rank]
        by_element.append((sorted(loc.primary for loc in locs), sorted(int(loc) for loc in locs)))
    return primary, complete, by_element


def _local_unsaturation(amp, numbered, exocyclic):
    """(compound locant count, hydro/ene/yne locants as a set, double bond locants) of one amplificant numbering
    (P-31.1.4.2, P-31.1.4.3)."""
    loc = numbered.local
    bonds = [sorted((loc[a], loc[b])) + [order] for a, b, order in amp.bonds]
    compound = sum(1 for low, high, _ in bonds if int(high) - int(low) != 10)
    positions = [int(low) for low, _, _ in bonds] + [int(h) for h in numbered.hydro] + [int(loc[a]) for a in exocyclic]
    return compound, sorted(positions), sorted(int(low) for low, _, order in bonds if order == 2)


def _multiple_bonds(mol, phane, position, picks):
    """[(low locant, high locant, order, bond index)] of every double or triple bond cited by an ending."""

    def end(node, atom):
        kind, value = phane.nodes[node]
        if kind == "amp":
            return PhaneLoc(position[node], picks[node].local[atom])
        return PhaneLoc(position[node])

    found = []
    for node, (kind, value) in enumerate(phane.nodes):
        if kind == "amp":
            for a, b, order in phane.amplificants[value].bonds:
                found.append((end(node, a), end(node, b), order, mol.GetBondBetweenAtoms(a, b).GetIdx()))
    for i, a, j, b, order, index in phane.multiple:
        found.append((end(i, a), end(j, b), order, index))
    return [(*sorted((x, y), key=int), order, index) for x, y, order, index in found]


def _is_simple(low, high):
    if low.local is None and high.local is None:
        return high.primary - low.primary == 1
    return low.local is not None and high.local is not None and low.primary == high.primary and int(high.local) - int(low.local) == 1000


def _evaluate(mol, phane, skeleton, free_atom, suffix_roots, stereo):
    name, outer, order = skeleton
    position = {node: i + 1 for i, node in enumerate(order)}
    free = _free_atoms(free_atom)
    seniority = {r: i for i, r in enumerate(sorted({a.rank for a in phane.amplificants}))}
    chain_position = {atom: position[k] for k, (kind, atom) in enumerate(phane.nodes) if kind == "atom"}
    ring_branches = {k: [] for k in phane.attachments}
    chain_branches = []
    owner_node = {}
    for k, (kind, value) in enumerate(phane.nodes):
        if kind == "amp":
            for atom in phane.amplificants[value].atoms:
                owner_node[atom] = k
    for owner, root, _ in phane.branches:
        if owner in owner_node:
            ring_branches[owner_node[owner]].append((owner, root))
        else:
            chain_branches.append((PhaneLoc(chain_position[owner]), owner, root))
    chosen, picks = {}, {}
    subs = list(chain_branches)
    free_locs = [PhaneLoc(chain_position[a]) for a in free if a in chain_position]
    added_loc = None
    ih_locs = []
    hetero = [
        (_RANK.get(mol.GetAtomWithIdx(atom).GetSymbol(), 99), PhaneLoc(chain_position[atom]))
        for atom in chain_position
        if mol.GetAtomWithIdx(atom).GetAtomicNum() != 6
    ]
    for k, (kind, value) in enumerate(phane.nodes):
        if kind != "amp":
            continue
        amp = phane.amplificants[value]
        cited = sorted(phane.attachments[k], key=lambda t: position[t[1]])
        owners = ring_branches[k]
        exocyclic = [a for i, a, j, b, _, _ in phane.multiple for node, a in ((i, a), (j, b)) if node == k]

        def local_key(numbered):
            loc = numbered.local
            return (
                sorted(loc[a] for a, _ in cited),
                [loc[a] for a, _ in cited],
                sorted(loc[a] for a in amp.replaced),
                sorted(numbered.ih),
                sorted(loc[o] for o, r in owners if r in suffix_roots),
                sorted(loc[a] for a in free if a in loc),
                [loc[amp.added]] if amp.added is not None else [],
                _local_unsaturation(amp, numbered, exocyclic),
                sorted(loc[o] for o, _ in owners),
            )

        numbered = min(amp.numberings, key=local_key)
        picks[k] = numbered
        loc = numbered.local
        p = position[k]
        chosen[k] = (sorted(loc[a] for a, _ in cited), [loc[a] for a, _ in cited])
        subs.extend((PhaneLoc(p, loc[o]), o, r) for o, r in owners)
        free_locs.extend(PhaneLoc(p, loc[a]) for a in free if a in loc)
        ih_locs.extend(PhaneLoc(p, h) for h in numbered.ih)
        if amp.added is not None:
            added_loc = PhaneLoc(p, loc[amp.added])
        hetero.extend((_RANK.get(symbol, 99), PhaneLoc(p, loc[a])) for a, symbol in amp.replaced.items())
    amp_nodes = sorted(k for k, node in enumerate(phane.nodes) if node[0] == "amp")
    superatoms = sorted(position[k] for k in amp_nodes)
    at = {position[k]: k for k in amp_nodes}
    ranks = {position[k]: seniority[phane.amplificants[phane.nodes[k][1]].rank] for k in amp_nodes}
    prefix_of = {position[k]: phane.amplificants[phane.nodes[k][1]].prefix for k in amp_nodes}
    singles = [p for p in superatoms if sum(1 for q in superatoms if prefix_of[q] == prefix_of[p]) == 1]
    by_seniority = sorted(set(prefix_of.values()), key=lambda s: min(ranks[p] for p in superatoms if prefix_of[p] == s))
    hetero_primary, hetero_complete, hetero_elements = _hetero_key(hetero)
    suffix_locs = sorted((loc for loc, _, root in subs if root in suffix_roots), key=int)
    bonds = _multiple_bonds(mol, phane, position, picks)
    hydro_locs = [PhaneLoc(position[k], h) for k, numbered in picks.items() for h in numbered.hydro]
    unsaturation = (
        sum(1 for low, high, _, _ in bonds if not _is_simple(low, high)),
        sorted([int(low) for low, _, _, _ in bonds] + [int(h) for h in hydro_locs]),
        sorted(int(low) for low, _, order, _ in bonds if order == 2),
    )
    descriptors = _descriptors(phane, position, picks, bonds, stereo)
    tiers = [[p for p in superatoms if ranks[p] == t] for t in sorted(set(ranks.values()))]
    citation = [p for s in by_seniority for p in superatoms if prefix_of[p] == s]
    # P-44.2.2.2.2 (c)-(h), P-44.2.2.2.6 (d)-(h): a linear phane compares the senior amplificants before all of them
    amplificant_key = (superatoms, tiers, citation) if phane.cyclic else (tiers, superatoms, citation)
    key = (
        outer,
        *amplificant_key,
        [[chosen[at[p]][0] for p in superatoms if prefix_of[p] == s] for s in by_seniority],
        [chosen[at[p]][1] for p in sorted(singles, key=lambda p: ranks[p])],
        [chosen[at[p]][1] for p in superatoms],
        hetero_primary,
        hetero_complete,
        hetero_elements,
        sorted(int(h) for h in ih_locs),
        [int(loc) for loc in suffix_locs],
        sorted((loc.primary, int(loc)) for loc in free_locs),
        [int(added_loc)] if added_loc is not None else [],
        unsaturation,
        sorted(loc.primary for loc, _, _ in subs),
        sorted(int(loc) for loc, _, _ in subs),
        [_CODE_RANK[d[2]] for d in descriptors],
    )
    return key, name, position, chosen, subs, free_locs, added_loc, ih_locs, hetero, bonds, hydro_locs, descriptors


def _descriptors(phane, position, picks, bonds, stereo):
    """[(locant, cited text, code, element)] of the specified stereo elements lying on the skeleton, the amplificants
    or an 'ene' bond, by ascending locant (P-31.1.4.3.4, P-91.3); the others belong to substituent groups."""
    where = {}
    for k, (kind, value) in enumerate(phane.nodes):
        if kind == "atom":
            where[value] = PhaneLoc(position[k])
        else:
            where.update((atom, PhaneLoc(position[k], local)) for atom, local in picks[k].local.items() if atom in phane.amplificants[value].atoms)
    cited = {index: str(low) if _is_simple(low, high) else f"{low}({high})" for low, high, order, index in bonds if order == 2}
    low_of = {index: low for low, _, order, index in bonds if order == 2}
    found = []
    for kind, index, code in stereo:
        if kind == "atom" and index in where:
            found.append((where[index], str(where[index]), code, (kind, index)))
        elif kind == "bond" and index in cited:
            found.append((low_of[index], cited[index], code, (kind, index)))
    return sorted(found, key=lambda t: int(t[0]))


def _amplification_text(phane, position, chosen):
    groups = {}
    for k, (kind, value) in enumerate(phane.nodes):
        if kind == "amp":
            groups.setdefault(phane.amplificants[value].prefix, []).append((position[k], tuple(chosen[k][1])))
    ranked = sorted(groups, key=lambda s: min(phane.amplificants[v].rank for kind, v in phane.nodes if kind == "amp" and phane.amplificants[v].prefix == s))
    pieces = []
    for prefix in ranked:
        patterns = {}
        for p, att in groups[prefix]:
            patterns.setdefault(att, []).append(p)
        text = ",".join(
            f"{','.join(map(str, sorted(locs)))}({','.join(map(str, att))})" for att, locs in sorted(patterns.items(), key=lambda kv: min(kv[1]))
        )
        amount = len(groups[prefix])
        if amount == 1:
            pieces.append(f"{text}-{prefix}")
        elif needs_bis(prefix):
            pieces.append(f"{text}-{multiplying_prefix(amount, True)}({prefix})")
        else:
            pieces.append(f"{text}-{multiplying_prefix(amount)}{prefix}")
    return "-".join(pieces)


def _replacement_text(hetero):
    by_element = {}
    for rank, loc in hetero:
        by_element.setdefault(rank, []).append(loc)
    symbol_of = {rank: symbol for symbol, rank in _RANK.items()}
    pieces = []
    for rank in sorted(by_element):
        locs = sorted(by_element[rank], key=int)
        word = _PREFIX.get(symbol_of.get(rank))
        if word is None:
            raise UnsupportedStructure("an unsupported replacement atom in a phane skeleton")
        pieces.append(f"{','.join(map(str, locs))}-{multiplying_prefix(len(locs)) if len(locs) > 1 else ''}{word}")
    return "-".join(pieces)


def phane_seniority_key(mol, phane=None):
    """P-44.2.2.2.2 (cyclic) or P-44.2.2.2.6 (linear) key of the phane `mol` (a smaller key is senior), or None when
    `mol` is not a parent phane of P-52.2.5.1."""
    phane = phane or find_phane(mol)
    if phane is None or not is_pin_phane(phane):
        return None
    ranked = [(skeleton[0], _evaluate(mol, phane, skeleton, None, frozenset(), [])[0]) for skeleton in phane.skeleton]
    name, best = min(ranked, key=lambda k: k[1])
    amplificants = tuple(sorted(a.rank for a in phane.amplificants))
    if phane.cyclic:
        skeleton_type = 0 if name.startswith(("spiro", "dispiro", "trispiro")) else 1 if "cyclo[" in name else 2
        return (skeleton_type, amplificants, best)
    hetero = [mol.GetAtomWithIdx(a).GetSymbol() for kind, a in phane.nodes if kind == "atom"]
    hetero = [symbol for symbol in hetero if symbol != "C"] + [s for a in phane.amplificants for s in a.replaced.values()]
    counts = tuple(-sum(1 for symbol in hetero if _RANK.get(symbol, 99) == r) for r in sorted(set(_RANK.values())))
    # (a) the senior amplificant, (b) most amplificants in order of seniority, (c) most nodes, (d)-(h) locants, (i)-(j)
    return (amplificants[:1], (*amplificants, ()), -phane.size, *best[:7], -len(hetero), counts, best[7:])


def name_phane_general(mol, free_atom=None, free_order=1, ignore=frozenset()):
    phane = find_phane(mol, free_atom, free_order, ignore)
    if phane is None:
        raise UnsupportedStructure("this structure is not a supported phane")
    stereo = specified_stereo_elements(mol) or []
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    suffix_class, suffix_roots = _principal_suffix(mol, phane.branches) if free_atom is None else (None, frozenset())
    valence = len(_free_atoms(free_atom))
    best = None
    for skeleton in phane.skeleton:
        result = _evaluate(mol, phane, skeleton, free_atom, suffix_roots, stereo)
        if best is None or result[0] < best[0]:
            best = result
    _, skeleton_name, position, chosen, subs, free_locs, added_loc, ih_locs, hetero, bonds, hydro_locs, descriptors = best
    suffix_locs = sorted((loc for loc, _, root in subs if root in suffix_roots), key=int)
    grouped = {}
    located = {d[3] for d in descriptors}
    context = {
        "atoms": {i: code for kind, i, code in stereo if kind == "atom" and (kind, i) not in located},
        "bonds": {
            (mol.GetBondWithIdx(i).GetBeginAtomIdx(), mol.GetBondWithIdx(i).GetEndAtomIdx()): code
            for kind, i, code in stereo
            if kind == "bond" and (kind, i) not in located
        },
        "used": set(),
    }
    token = BRANCH_STEREO.set(context)
    try:
        for loc, owner, root in subs:
            if root in suffix_roots:
                continue
            name, compound = name_branch(graph, root, owner, halogens, frozenset(), mol=mol, unsaturated=True)
            grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(loc)
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element of a phane substituent is not cited by any supported name")
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    hydro = (
        f"{','.join(map(str, sorted(hydro_locs, key=int)))}-{multiplied_word(len(hydro_locs), 'hydro')}" if hydro_locs else ""
    )
    amplification = _amplification_text(phane, position, chosen)
    replacement = _replacement_text(hetero) if hetero else ""
    indicated = ",".join(f"{loc}H" for loc in sorted(ih_locs, key=int)) + "-" if ih_locs else ""
    core = indicated + (replacement + "-" if replacement else "") + amplification + skeleton_name
    ordered = sorted(bonds, key=lambda b: (int(b[0]), int(b[1])))
    cited = {index: str(low) if _is_simple(low, high) else f"{low}({high})" for low, high, _, index in ordered}
    segments = []
    for order, word in ((2, "ene"), (3, "yne")):
        locants = [cited[index] for _, _, o, index in ordered if o == order]
        if locants:
            segments.append((",".join(locants), multiplied_word(len(locants), word)))
    if free_atom is not None:
        added_text = f"({added_loc}H)" if added_loc is not None else ""
        word = "ylidene" if free_order == 2 else multiplied_word(valence, "yl")
        segments.append((f"{','.join(str(loc) for loc in sorted(free_locs, key=int))}{added_text}", word))
    if suffix_locs:
        segments.append((",".join(str(loc) for loc in suffix_locs), multiplied_word(len(suffix_locs), _SUFFIX_WORDS[suffix_class])))
    for locants, word in segments:
        core = (core[:-1] if core.endswith("e") and elides_before(word) else core) + f"-{locants}-{word}"
    name = "-".join(part for part in (prefix, hydro, core) if part)
    return ("(" + ",".join(f"{d[1]}{d[2]}" for d in descriptors) + ")-" if descriptors else "") + name


def phane_substituent(mol, graph, root, coming_from, preferred=True):
    """(name, True) of a phane group entered at `root`, or None when the branch is not a supported phane; with
    `preferred` False, the non-preferred phane group of a fused cyclophane."""
    arm = _component(graph, root, {coming_from})
    order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
    if order not in (1.0, 2.0) or sum(1 for ring in mol.GetRingInfo().AtomRings() if set(ring) <= arm) < 2:
        return None
    editable = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - arm, reverse=True):
        editable.RemoveAtom(index)
    kept = sorted(arm)
    free_atom = kept.index(root)
    ignore = frozenset()
    if order == 2.0:
        dummy = editable.AddAtom(Chem.Atom(0))
        editable.AddBond(free_atom, dummy, Chem.BondType.DOUBLE)
        ignore = frozenset({dummy})
    fragment = editable.GetMol()
    try:
        sanitize_probe(fragment)
    except Exception:
        return None
    found = find_phane(fragment, free_atom, int(order), ignore)
    if found is None or not (is_pin_phane(found) if preferred else is_fused_cyclophane(found)):
        return None
    name = name_phane_general(fragment, free_atom, int(order), ignore)
    if not preferred:
        mark(name, "P-52.2.5.2 prefers a fusion or bridged fused group to this phane group")
    return name, True


def phane_diyl_name(mol, skeleton_atoms, free_atoms, blocked):
    """Name of the phane made of `skeleton_atoms` (plus the substituents hanging on it) with a free valence at each of
    `free_atoms` (P-29.3.6, P-29.2); `blocked` are the atoms of the units joined through those valences. None when
    the skeleton is not a supported phane."""
    keep = set(skeleton_atoms)
    graph = adjacency(mol)
    stack = list(keep)
    while stack:
        for n in graph[stack.pop()]:
            if n not in keep and n not in blocked:
                keep.add(n)
                stack.append(n)
    editable = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(index)
    kept = sorted(keep)
    fragment = editable.GetMol()
    try:
        Chem.SanitizeMol(fragment)
    except Exception:
        return None
    free = tuple(kept.index(a) for a in free_atoms)
    if find_phane(fragment, free, 1) is None:
        return None
    return name_phane_general(fragment, free, 1)
