"""P-25.3 fusion names of ortho- and ortho- and peri-fused mancude ring systems: the ring system is divided into a
parent component and first- and higher-order attached components (P-25.3.2.4, P-25.3.4), the fusion descriptors are
built from the peripheral lettering of the parent and the locants of the attached components (P-25.3.1.3, P-25.3.4.2.4,
P-25.3.6, P-25.3.8) and the whole system is numbered by P-25.3.3."""

import re
from contextvars import ContextVar
from dataclasses import dataclass, field
from itertools import combinations, permutations, product

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure
from ._fused_numbering import FusedSystem, _locant_key, fused_numberings
from ._fusion_components import Component, exception_numbering, identify, skeleton
from ._numerals import numerical_term
from ._pin import mark

FUSION_NAME_REQUIRED = ContextVar("fusion_name_required", default=False)
PREFER_VON_BAEYER = ContextVar("prefer_von_baeyer", default=False)

_MAX_RINGS_IN_COMPONENT = 10
_MAX_RINGS_IN_SYSTEM = 14


@dataclass
class Part:
    rings: frozenset
    atoms: list
    comp: Component
    sub: Chem.Mol
    parent: "Part" = None
    children: list = field(default_factory=list)
    order: int = 0
    role: str = "attached"

    def numberings_in_system(self):
        return [{self.atoms[i]: loc for i, loc in n.items()} for n in self.comp.numberings]


class Context:
    def __init__(self, mol):
        Chem.GetSymmSSSR(mol)
        self.mol = mol
        self.sk = skeleton([a.GetSymbol() for a in mol.GetAtoms()], [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])
        self.fs = FusedSystem(self.sk)
        self.n = len(self.fs.rings)
        self.adj = {i: set(self.fs.adj[i]) for i in range(self.n)}
        self._parts = {}
        self.subsets = self._connected_subsets()

    def _connected_subsets(self):
        found = set()
        frontier = {frozenset([i]) for i in range(self.n)}
        found |= frontier
        while frontier:
            nxt = set()
            for s in frontier:
                if len(s) >= _MAX_RINGS_IN_COMPONENT:
                    continue
                for r in s:
                    for nb in self.adj[r]:
                        if nb not in s:
                            t = s | {nb}
                            if t not in found:
                                found.add(t)
                                nxt.add(t)
            frontier = nxt
        return sorted(found, key=lambda s: (len(s), sorted(s)))

    def parts_for(self, subset):
        if subset in self._parts:
            return self._parts[subset]
        atoms = sorted(set().union(*[self.fs.ring_sets[r] for r in subset]))
        index = {a: i for i, a in enumerate(atoms)}
        bonds = set()
        for r in subset:
            cyc = self.fs.rings[r]
            for k in range(len(cyc)):
                bonds.add(frozenset((cyc[k], cyc[(k + 1) % len(cyc)])))
        sub = skeleton([self.sk.GetAtomWithIdx(a).GetSymbol() for a in atoms], [tuple(index[x] for x in b) for b in bonds])
        parts = []
        if sub.GetRingInfo().NumRings() == len(subset):
            try:
                comps = identify(sub)
            except UnsupportedStructure:
                comps = []
            parts = [Part(subset, atoms, c, sub) for c in comps]
        self._parts[subset] = parts
        return parts

    def shared(self, a_rings, b_rings):
        return any(j in b_rings for i in a_rings for j in self.adj[i])


def _tree_decompositions(ctx, roots, parents=(), initial=None):
    """Every division of the ring system into valid components forming a tree. `roots` is the parent component, or the
    interparent component of a multiparent system whose `parents` (parts) are attached to it; `initial` places the parts
    of an extended multiparent system as (part, index of its parent part, role)."""
    results = []
    start = initial or [(roots, None, "inter" if parents else "root")] + [(p, 0, "parent") for p in parents]

    def build(placed):
        parts = []
        for proto, parent_index, role in placed:
            parent = parts[parent_index] if parent_index is not None else None
            base, _, explicit = role.partition(":")
            order = int(explicit) if explicit else {"parent": 0, "inter": 1, "center": 2}.get(base, 0 if not parent else parent.order + 1)
            part = Part(proto.rings, proto.atoms, proto.comp, proto.sub, parent, [], order, base)
            if parent:
                parent.children.append(part)
            parts.append(part)
        return parts

    def recurse(placed, unplaced):
        if not unplaced:
            results.append(build(placed))
            return
        placed_rings = set().union(*[p.rings for p, _, _ in placed])
        candidates = [r for r in unplaced if any(nb in placed_rings for nb in ctx.adj[r])]
        frontier = min(candidates)
        for subset in ctx.subsets:
            if frontier not in subset or not subset <= unplaced:
                continue
            for proto in ctx.parts_for(subset):
                touching = [i for i, (p, _, _) in enumerate(placed) if ctx.shared(p.rings, subset)]
                if len(touching) != 1:
                    continue
                placed.append((proto, touching[0], "attached"))
                recurse(placed, unplaced - subset)
                placed.pop()

    used = set().union(*[p.rings for p, _, _ in start])
    recurse(start, frozenset(range(ctx.n)) - used)
    return results


def _signature(part):
    """Identity of a part with its subtree, for recognising identical attached components."""
    return (part.comp.name, part.comp.hetero_text, tuple(sorted(_signature(c) for c in part.children)))


# P-25.3.5.3: a multiparent name beats the retained benzoazole names (benzo[1,2-c:4,5-c']dipyrrole, P-72.2, P-73.5.1.3)
_BENZOAZOLES = {"indole", "isoindole", "indazole"}


def _valid_decomposition(ctx, parts, retained_sets):
    """P-25.3.5: retained components are never broken up by weaker ones, and an isolated benzene ring fused to a
    heteromonocycle of five or more members belongs to a benzoheterocycle unit, unless that would break up a multiplicative
    prefix or a multiparent name."""
    fs, sk = ctx.fs, ctx.sk
    part_of = {r: p for p in parts for r in p.rings}
    weak = lambda p: p.comp.kind == "mono" or p.comp.benzo_unit
    for ring_set in retained_sets:
        owners = {id(part_of[r]): part_of[r] for r in ring_set}
        if len(owners) > 1 and all(weak(p) or p.rings <= ring_set for p in owners.values()):
            return False
    elements = lambda r: {sk.GetAtomWithIdx(a).GetSymbol() for a in fs.rings[r]}
    isolated = lambda r: len(fs.rings[r]) == 6 and elements(r) == {"C"} and weak(part_of[r])
    benzene = {r for r in range(ctx.n) if isolated(r)}
    hosts = {r: {b for b in ctx.adj[r] if b in benzene} for r in range(ctx.n) if len(fs.rings[r]) >= 5 and elements(r) != {"C"}}
    kind_of = lambda h: tuple(sorted(sk.GetAtomWithIdx(a).GetSymbol() for a in fs.rings[h]))
    host_count = {}
    for r in benzene:
        kinds = [kind_of(h) for h, bs in hosts.items() if r in bs]
        host_count[r] = max((kinds.count(k) for k in set(kinds)), default=0)
    for part in parts:
        comp = part.comp
        if comp.benzo_unit:
            rings = sorted(part.rings)
            benzo = [r for r in rings if elements(r) == {"C"}]
            host = [r for r in rings if r in hosts]
            if not benzo or not host or benzo[0] not in benzene:
                return False
            if len(hosts[host[0]]) > 1 or host_count.get(benzo[0], 0) > 1:
                return False
        elif comp.kind == "mono" and len(part.rings) == 1:
            (r,) = part.rings
            if r in benzene and host_count[r] < 2:
                for h, bs in hosts.items():
                    if r in bs and len(bs) == 1:
                        owner = part_of[h]
                        if owner.comp.kind == "mono" and owner is not part:
                            return False
    return True


def _seniority_key(parts):
    """P-25.3.4.2.2 (a): the senior ring system at first order, then at second order, and so on."""
    depth = max(p.order for p in parts)
    return tuple(tuple(sorted(p.comp.senior_key for p in parts if p.order == k and p.role != "parent")) for k in range(1, depth + 1))


def _location_key(parts):
    """P-25.3.4.2.1 (b)-(d): fewest orders, most attached components per order, most identical components multiplied."""
    depth = max(p.order for p in parts)
    counts = [sum(1 for p in parts if p.order == k and p.role != "parent") for k in range(1, depth + 1)]
    multiplied = 0
    for p in parts:
        groups = {}
        for c in p.children:
            if c.role != "parent":
                groups.setdefault(_signature(c), []).append(c)
        multiplied += sum(len(g) for g in groups.values() if len(g) > 1)
    return (depth, tuple(-c for c in counts), -multiplied)


def _location_of(parts):
    return frozenset(p.rings for p in parts if p.role in ("root", "parent", "inter"))


def _preferred(decomps):
    """The decompositions left after choosing the location of the parent (P-25.3.4.2.1) and then the attached components
    (P-25.3.4.2.2 (a))."""
    by_location = {}
    for d in decomps:
        by_location.setdefault(_location_of(d), []).append(d)
    best_per_location = []
    for group in by_location.values():
        top = min(_seniority_key(d) for d in group)
        best_per_location.extend(d for d in group if _seniority_key(d) == top)
    top = min(_location_key(d) for d in best_per_location)
    return [d for d in best_per_location if _location_key(d) == top]


def _periphery_sequence(part, numbering_sys):
    """Peripheral atoms of the part in the direction of its numbering, starting at locant 1."""
    sub_fs = FusedSystem(part.sub)
    cycle = [part.atoms[i] for i in sub_fs.periphery]
    locant_of = numbering_sys
    first = next(a for a in cycle if locant_of[a] == "1")
    pos = cycle.index(first)
    forward = cycle[pos:] + cycle[:pos]
    backward = [forward[0]] + forward[:0:-1]
    nxt_f, nxt_b = locant_of[forward[1]], locant_of[backward[1]]
    return forward if _locant_key(nxt_f) <= _locant_key(nxt_b) else backward


def _letter(index):
    letters = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        letters = chr(97 + rem) + letters
    return letters


def _path_on_periphery(sequence, shared_atoms):
    """The shared atoms as a contiguous stretch of the cyclic periphery `sequence`, in order, or None."""
    n = len(sequence)
    flags = [a in shared_atoms for a in sequence]
    if all(flags):
        return None
    for start in range(n):
        if flags[start] and not flags[start - 1]:
            path = []
            k = start
            while flags[k % n] and len(path) < n:
                path.append(sequence[k % n])
                k += 1
            if len(path) == len(shared_atoms):
                return path, start
    return None


def _prime(locant, count):
    return f"{locant}{chr(39) * count}"


@dataclass
class Rendered:
    text: str
    key: tuple


def _sources_items(sources):
    """Attached components of the given (part, numbering, letter primes) sources with their stretch of the periphery."""
    items = []
    for part, numbering_sys, letter_primes in sources:
        if not part.children:
            continue
        sequence = _periphery_sequence(part, numbering_sys)
        letters_of = {frozenset((sequence[k], sequence[(k + 1) % len(sequence)])): k for k in range(len(sequence))}
        for child in part.children:
            if child.role == "parent":
                continue
            shared = set(child.atoms) & set(part.atoms)
            located = _path_on_periphery(sequence, shared)
            if located is None:
                return None
            path, _ = located
            bond_indices = [letters_of[frozenset((path[i], path[i + 1]))] for i in range(len(path) - 1)]
            items.append((child, path, bond_indices, letter_primes, part, numbering_sys))
    return items


def _render_children(ctx, sources, parent_primes, is_root):
    """Prefix text and rank key for the attached components of the (part, numbering, letter primes) sources."""
    items = _sources_items(sources)
    if items is None:
        return None
    if not items:
        return Rendered("", ())
    groups = {}
    for item in items:
        groups.setdefault(_signature(item[0]), []).append(item)
    rendered_groups = []
    for members in groups.values():
        rendered = _render_group(ctx, members, parent_primes, is_root)
        if rendered is None:
            return None
        rendered_groups.extend(rendered)
    rendered_groups.sort(key=lambda r: (r["name"], r["hetero"], r["order_key"]))
    letters = tuple(sorted(i for r in rendered_groups for member in r["rank"][0] for i in member))
    return Rendered(_join_prefixes([r["text"] for r in rendered_groups]), (letters,) + tuple(r["rank"] for r in rendered_groups))


def _join_prefixes(texts):
    out = ""
    for text in texts:
        if out and re.match(r"(?:as|s)-", text):
            out += "-"
        out += text
    return out


def _locant_cited(child, path, child_numbering):
    locants = [child_numbering[a] for a in path]
    if len(locants) > 2:
        locants = [loc for loc in locants if loc.isdigit()]
    return locants


def _prime_letters(letters, count):
    return "".join(c + chr(39) * count for c in letters) if count else letters


def _render_group(ctx, members, parent_primes, is_root):
    first = members[0][0]
    source_part = members[0][4]
    child_level = first.order
    hydrocarbon_monocycle = first.comp.kind in ("mono", "attached_only") and all(e == "C" for e in first.comp.elements)
    omit_own = hydrocarbon_monocycle and not first.children
    bare = (
        is_root and omit_own and len(members) == 1 and source_part.comp.kind == "mono" and all(e == "C" for e in source_part.comp.elements)
        and source_part.role == "root"
    )
    best_per_member = []
    for child, path, bond_indices, letter_primes, src_part, src_numbering in members:
        options = []
        for numbering in child.numberings_in_system():
            cited = _locant_cited(child, path, numbering)
            sub = _render_children(ctx, [(child, numbering, 0)], 0, False)
            if sub is None:
                continue
            options.append((tuple(sorted(_locant_key(l) for l in cited)), tuple(_locant_key(l) for l in cited), sub.key, cited, sub, numbering))
        if not options:
            return None
        best_per_member.append((child, path, bond_indices, min(options, key=lambda o: o[:3]), letter_primes, src_part, src_numbering))
    best_per_member.sort(key=lambda m: (m[2], m[4], m[3][0], m[3][1]))
    if len(best_per_member) > 1 and first.children:
        nested_texts = {_render_children(ctx, [(m[0], m[3][5], 0)], m[0].order - 1, False).text for m in best_per_member}
        if len(nested_texts) > 1:
            singles = []
            for member in members:
                rendered = _render_group(ctx, [member], parent_primes, is_root)
                if rendered is None:
                    return None
                singles.extend(rendered)
            return singles
    descriptors, resolved = [], []
    for index, (child, path, bond_indices, best, letter_primes, src_part, src_numbering) in enumerate(best_per_member):
        cited = best[3]
        own_primes = (child_level - 1) + index
        resolved.append(_render_children(ctx, [(child, best[5], 0)], own_primes, False))
        if is_root:
            letters = _prime_letters("".join(_letter(i) for i in bond_indices), letter_primes)
            if omit_own:
                descriptors.append((letters, "abbrev"))
            else:
                descriptors.append((",".join(_prime(l, own_primes) for l in cited) + "-" + letters, "complete"))
        else:
            parent_loc = [src_numbering[a] for a in path]
            if len(path) > 2:
                parent_loc = [l for l in parent_loc if l.isdigit()]
            plocs = ",".join(_prime(l, parent_primes) for l in parent_loc)
            if omit_own:
                descriptors.append((plocs, "abbrev"))
            else:
                descriptors.append((",".join(_prime(l, own_primes) for l in cited) + ":" + plocs, "complete"))
    multiplied = len(best_per_member) > 1
    if multiplied:
        complete = descriptors[0][1] == "complete"
        separator = (":" if is_root else ";") if complete else (":" if not is_root else ",")
        joined = separator.join(d[0] for d in descriptors)
    else:
        joined = descriptors[0][0]
    prefix_core = first.comp.bracketed_prefix()
    subtext = resolved[0].text
    if multiplied:
        count = len(best_per_member)
        mult = numerical_term(count)
        if first.children or prefix_core.startswith("["):
            mult = {2: "bis", 3: "tris", 4: "tetrakis"}.get(count, mult)
            text = mult + "(" + subtext + prefix_core + ")" + f"[{joined}]"
        else:
            text = subtext + mult + prefix_core + f"[{joined}]"
    elif bare:
        text = subtext + prefix_core
    else:
        text = subtext + prefix_core + f"[{joined}]"
    rank = tuple((m[3][0], m[3][1], m[3][2]) for m in best_per_member)
    citation = re.sub(r"\[[^\]]*\]|[()]", "", subtext + prefix_core).lower()
    return [
        {
            "text": text,
            "name": citation,
            "hetero": first.comp.hetero_text,
            "order_key": rank[0],
            "rank": (tuple(tuple(m[2]) for m in best_per_member), tuple(m[4] for m in best_per_member), rank),
        }
    ]


def _multiparent_decompositions(ctx, group):
    """P-25.3.7.1: two or more nonoverlapping locations of the parent component fused to one first-order interparent
    component."""
    by_type = {}
    for part in group:
        by_type.setdefault((part.comp.name, part.comp.hetero_text), []).append(part)
    results = []
    for placements in by_type.values():
        for m in range(2, len(placements) + 1):
            for combo in combinations(placements, m):
                if any(a.rings & b.rings or ctx.shared(a.rings, b.rings) for a, b in combinations(combo, 2)):
                    continue
                taken = set().union(*[p.rings for p in combo])
                for subset in ctx.subsets:
                    if subset & taken or not all(ctx.shared(subset, p.rings) for p in combo):
                        continue
                    for inter in ctx.parts_for(subset):
                        results.extend(_tree_decompositions(ctx, inter, combo))
    return results


def _render_multiparent(ctx, parts):
    """Name text in front of the multiplied parent, and its rank key, for a multiparent decomposition."""
    inter = parts[0]
    parents = [p for p in parts if p.role == "parent"]
    letters_needed = not (
        parents[0].comp.kind == "mono" and all(e == "C" for e in parents[0].comp.elements) and not any(p.children for p in parents)
    )
    best = None
    for ordering in set(permutations(range(len(parents)))):
        ordered = [parents[i] for i in ordering]
        for combo in product(*[p.numberings_in_system() for p in ordered]):
            for inter_numbering in inter.numberings_in_system():
                descriptors, keys = [], []
                valid = True
                for index, (parent, numbering) in enumerate(zip(ordered, combo)):
                    sequence = _periphery_sequence(parent, numbering)
                    located = _path_on_periphery(sequence, set(inter.atoms) & set(parent.atoms))
                    if located is None:
                        valid = False
                        break
                    path, _ = located
                    letters_of = {frozenset((sequence[k], sequence[(k + 1) % len(sequence)])): k for k in range(len(sequence))}
                    indices = [letters_of[frozenset((path[i], path[i + 1]))] for i in range(len(path) - 1)]
                    cited = _locant_cited(inter, path, inter_numbering)
                    letters = _prime_letters("".join(_letter(i) for i in indices), index)
                    descriptors.append(",".join(cited) + ("-" + letters if letters_needed else ""))
                    keys.append((tuple(indices), tuple(_locant_key(l) for l in cited)))
                if not valid:
                    continue
                letter_key = (tuple(sorted(k[0] for k in keys)), tuple(k[0] for k in keys))
                locant_key = (tuple(sorted(x for k in keys for x in k[1])), tuple(k[1] for k in keys))
                sources = [(p, n, i) for i, (p, n) in enumerate(zip(ordered, combo))]
                on_parents = _render_children(ctx, sources, 0, True)
                on_inter = _render_children(ctx, [(inter, inter_numbering, 0)], 0, False)
                if on_parents is None or on_inter is None:
                    continue
                key = (letter_key, locant_key, on_parents.key, on_inter.key)
                if best is None or key < best[0]:
                    best = (key, on_parents.text + on_inter.text + inter.comp.bracketed_prefix() + "[" + ":".join(descriptors) + "]")
    return best


def _branch_chains(ctx, branch, center_rings, group_types, valid_subsets):
    """Ways to divide the rings of one branch into a chain of interparent components (the one next to the center first)
    ending in a parent component: lists of (subset, part) with the parent last."""
    found = []

    def extend(remaining, chain):
        previous = chain[-1][0] if chain else center_rings
        others = [s for s, _ in chain[:-1]]
        if chain and remaining in valid_subsets and ctx.shared(remaining, previous) and not ctx.shared(remaining, center_rings) and not any(ctx.shared(remaining, s) for s in others):
            for part in ctx.parts_for(remaining):
                if (part.comp.name, part.comp.hetero_text) in group_types:
                    found.append(chain + [(remaining, part)])
        for subset in valid_subsets:
            if subset == remaining or not subset <= remaining or not ctx.shared(subset, previous):
                continue
            if chain and (ctx.shared(subset, center_rings) or any(ctx.shared(subset, s) for s in others)):
                continue
            rest = remaining - subset
            if not _connected(ctx, rest):
                continue
            for part in ctx.parts_for(subset):
                extend(rest, chain + [(subset, part)])

    extend(branch, [])
    return found


def _connected(ctx, rings):
    rings = set(rings)
    if not rings:
        return False
    seen, frontier = set(), {next(iter(rings))}
    while frontier:
        seen |= frontier
        frontier = {n for r in frontier for n in ctx.adj[r] if n in rings} - seen
    return seen == rings


def _extended_decompositions(ctx, group):
    """P-25.3.7.3: parent components separated by an odd number of symmetrically arranged interparent components: a
    center joined to two or more identical chains of first-, second-, ... order interparent components ending in a parent."""
    group_types = {(part.comp.name, part.comp.hetero_text) for part in group}
    everything = frozenset(range(ctx.n))
    valid_subsets = {subset for subset in ctx.subsets if ctx.parts_for(subset)}
    results = []
    for center_rings in valid_subsets:
        rest = everything - center_rings
        if not rest or len(center_rings) == ctx.n:
            continue
        branches, seen = [], set()
        for ring in sorted(rest):
            if ring in seen:
                continue
            component, frontier = set(), {ring}
            while frontier:
                component |= frontier
                frontier = {n for r in frontier for n in ctx.adj[r] if n in rest} - component
            seen |= component
            branches.append(frozenset(component))
        if len(branches) < 2 or any(not ctx.shared(b, center_rings) for b in branches):
            continue
        per_branch = [_branch_chains(ctx, b, center_rings, group_types, valid_subsets) for b in branches]
        if any(not chains for chains in per_branch):
            continue
        for center in ctx.parts_for(center_rings):
            for combo in product(*per_branch):
                signature = [tuple((p.comp.name, p.comp.hetero_text) for _, p in chain) for chain in combo]
                if any(sig != signature[0] for sig in signature) or len(signature[0]) < 2:
                    continue
                k = len(signature[0]) - 1
                initial = [(center, None, f"center:{k + 1}")]
                for chain in combo:
                    link = 0
                    for level, (subset, part) in zip(range(k, 0, -1), chain[:-1]):
                        initial.append((part, link, f"inter:{level}"))
                        link = len(initial) - 1
                    initial.append((chain[-1][1], link, "parent"))
                results.extend(_tree_decompositions(ctx, None, initial=initial))
    return results


def _chain_of(parent):
    """The interparent components from a parent up to the center, nearest the parent first."""
    chain = []
    part = parent.parent
    while part is not None and part.role == "inter":
        chain.append(part)
        part = part.parent
    return chain, part


# P-25.3.7.3 (a) gives primes by order only; chains with k >= 2 have no Blue Book example.
def _prime_level(level, branch):
    return branch if level == 1 else level


def _branch_descriptor(piece, numberings, center, n_center):
    """Letters, cited locants and, for every higher level, the own and target locants of one branch (parent first, then
    its interparent components up to the one next to the center) in the given numberings, or None."""
    parent, level_one = piece[0], piece[1]
    sequence = _periphery_sequence(parent, numberings[0])
    located = _path_on_periphery(sequence, set(level_one.atoms) & set(parent.atoms))
    if located is None:
        return None
    path, _ = located
    letters_of = {frozenset((sequence[i], sequence[(i + 1) % len(sequence)])): i for i in range(len(sequence))}
    indices = tuple(letters_of[frozenset((path[i], path[i + 1]))] for i in range(len(path) - 1))
    cited = _locant_cited(level_one, path, numberings[1])
    upper = []
    for level in range(2, len(piece) + 1):
        lower, lower_numbering = piece[level - 1], numberings[level - 1]
        above, above_numbering = (piece[level], numberings[level]) if level < len(piece) else (center, n_center)
        sequence_l = _periphery_sequence(lower, lower_numbering)
        located_l = _path_on_periphery(sequence_l, set(above.atoms) & set(lower.atoms))
        if located_l is None:
            return None
        path_l, _ = located_l
        targets = [lower_numbering[a] for a in path_l]
        owns = [above_numbering[a] for a in path_l]
        if len(path_l) > 2:
            targets = [t for t in targets if t.isdigit()]
            owns = [o for o in owns if o.isdigit()]
        upper.append((owns, targets))
    key = (indices, tuple(_locant_key(l) for l in cited), tuple(tuple(_locant_key(o) for o in owns) for owns, _ in upper))
    return key, indices, cited, upper


def _render_extended(ctx, parts):
    center = parts[0]
    parents = [p for p in parts if p.role == "parent"]
    chains = [_chain_of(p)[0] for p in parents]
    if len(parents) < 2 or any(not c for c in chains) or len({len(c) for c in chains}) != 1:
        return None
    k = len(chains[0])
    pieces = [(parent, *chain) for parent, chain in zip(parents, chains)]
    best = None
    for n_center in center.numberings_in_system():
        branches = []
        for piece in pieces:
            options = [
                described
                for numberings in product(*[part.numberings_in_system() for part in piece])
                if (described := _branch_descriptor(piece, numberings, center, n_center)) is not None
            ]
            if not options:
                branches = None
                break
            branches.append(min(options, key=lambda o: o[0]))
        if branches is None:
            continue
        branches.sort(key=lambda b: b[0])
        keys = [b[0] for b in branches]
        key = (
            (tuple(sorted(x[0] for x in keys)), tuple(x[0] for x in keys)),
            (tuple(sorted(y for x in keys for y in x[1])), tuple(x[1] for x in keys)),
            tuple(x[2] for x in keys),
        )
        if best is not None and key >= best[0]:
            continue
        inner, upper = [], [[] for _ in range(k)]
        for b, (_, indices, cited, levels) in enumerate(branches):
            inner.append(",".join(_prime(l, b) for l in cited) + "-" + _prime_letters("".join(_letter(i) for i in indices), b))
            for level, (owns, targets) in enumerate(levels, start=2):
                upper[level - 2].append(",".join(_prime(o, level) for o in owns) + ":" + ",".join(_prime(t, _prime_level(level - 1, b)) for t in targets))
        count = numerical_term(len(parents))
        text = center.comp.bracketed_prefix() + "[" + ";".join(upper[k - 1]) + "]"
        for level in range(k, 1, -1):
            text += count + pieces[0][level].comp.bracketed_prefix() + "[" + ";".join(upper[level - 2]) + "]"
        text += count + pieces[0][1].comp.bracketed_prefix() + "[" + ":".join(inner) + "]"
        best = (key, text, pieces[0][0], len(parents))
    return best


def _multiparent_name(parent_part, count):
    comp = parent_part.comp
    mult = numerical_term(count)
    if comp.hetero_text:
        return {2: "bis", 3: "tris", 4: "tetrakis"}.get(count, mult) + "(" + f"[{comp.hetero_text}]{comp.name}" + ")"
    return mult + comp.name


def _root_candidates(ctx, root_part):
    best = None
    for numbering in root_part.numberings_in_system():
        rendered = _render_children(ctx, [(root_part, numbering, 0)], 0, True)
        if rendered is None:
            continue
        key = rendered.key
        if best is None or key < best[0]:
            best = (key, rendered, numbering)
    return best


def _parent_name(root, attached):
    comp = root.comp
    if comp.hetero_text and attached:
        return f"[{comp.hetero_text}]{comp.name}"
    if comp.hetero_text:
        return f"{comp.hetero_text}-{comp.name}"
    return comp.name


def fusion_name(mol):
    """Fusion name (without indicated hydrogen) of the ortho- and peri-fused ring system `mol`, or UnsupportedStructure.
    A system with a third component ortho- and peri-fused to two others takes the skeletal replacement name of P-25.5.1."""
    name, root, _ = fusion_name_keyed(mol)
    return name, root


def fusion_name_keyed(mol):
    """`fusion_name` with the descriptor key of `_fusion_name_keyed` as a third value."""
    try:
        name, root, key = _fusion_name_keyed(mol)
    except UnsupportedStructure:
        if all(a.GetSymbol() == "C" for a in mol.GetAtoms()):
            raise
        name, root, key = _replacement_name(mol)
    if sum(len(ring) >= 5 for ring in Chem.GetSymmSSSR(mol)) < 2:
        mark(name, "fusion nomenclature gives preferred names only to systems with two rings of five or more members (P-52.2.4.1)")
    return name, root, key


def _replacement_name(mol):
    """P-25.5.1.1: the heteroatoms replace carbon atoms of the hydrocarbon fusion name; its numbering is unchanged."""
    from ._fusion_components import A_PREFIX
    from ._fused_numbering import HETERO_RANK

    skeleton_c = skeleton(["C"] * mol.GetNumAtoms(), [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])
    name, root, descriptor = _fusion_name_keyed(skeleton_c)
    options = system_numbering_options(Context(skeleton_c), name, root)
    hetero = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "C"]

    def key(numbering):
        by_locant = sorted(hetero, key=lambda h: _locant_key(numbering[h]))
        return (
            [_locant_key(numbering[h]) for h in by_locant],
            [HETERO_RANK[mol.GetAtomWithIdx(h).GetSymbol()] for h in by_locant],
        )

    best = min(key(n) for n in options)
    tied = [n for n in options if key(n) == best]
    numbering = tied[0]
    groups = {}
    for h in hetero:
        groups.setdefault(mol.GetAtomWithIdx(h).GetSymbol(), []).append(numbering[h])
    pieces = []
    for element in sorted(groups, key=lambda e: HETERO_RANK[e]):
        locants = sorted(groups[element], key=_locant_key)
        multiplier = "" if len(locants) == 1 else numerical_term(len(locants))
        pieces.append(",".join(locants) + "-" + multiplier + A_PREFIX[element])
    return "-".join(pieces) + name, _ReplacementRoot(root, tied), descriptor


class _ReplacementRoot:
    def __init__(self, root, numberings):
        self.root = root
        self.numberings = numberings
        self.rings = root.rings


def _fusion_name_core(mol):
    name, root, _ = _fusion_name_keyed(mol)
    return name, root


def _fusion_name_keyed(mol):
    """(name, root, descriptor key): the key holds the fusion descriptor letters and then the locants in order of
    appearance, `()` for a retained or multiparent name (P-44.2.2.2.3 (c)-(d))."""
    Chem.GetSymmSSSR(mol)
    ctx = Context(mol)
    whole = frozenset(range(ctx.n))
    for part in ctx.parts_for(whole):
        if part.comp.kind in ("hydro", "hetero"):
            return _parent_name(part, False), part, ()
    retained_sets = [
        subset
        for subset in ctx.subsets
        if len(subset) > 1 and any(p.comp.retained and p.comp.kind != "mono" and not p.comp.benzo_unit for p in ctx.parts_for(subset))
    ]
    roots = []
    for subset in ctx.subsets:
        for part in ctx.parts_for(subset):
            if part.comp.kind in ("mono", "hydro", "hetero"):
                roots.append(part)
    roots.sort(key=lambda p: p.comp.senior_key)
    index = 0
    while index < len(roots):
        group = [r for r in roots[index:] if r.comp.senior_key == roots[index].comp.senior_key]
        index += len(group)
        extended = [d for d in _extended_decompositions(ctx, group) if _valid_decomposition(ctx, d, retained_sets)]
        scored_ext = []
        for parts in extended:
            result = _render_extended(ctx, parts)
            if result is not None:
                scored_ext.append((result[0], result[1], parts, result[2], result[3]))
        if scored_ext:
            scored_ext.sort(key=lambda r: (r[0], r[1]))
            _, text, parts, parent0, count = scored_ext[0]
            return text + _multiparent_name(parent0, count), parts[0], ()
        multi_sets = [
            subset
            for subset in retained_sets
            if any(p.comp.retained and p.comp.kind != "mono" and not p.comp.benzo_unit and p.comp.name not in _BENZOAZOLES for p in ctx.parts_for(subset))
        ]
        multi_groups = [group]
        if group[0].comp.name in _BENZOAZOLES:
            later = roots[index:]
            keys = {repr(r.comp.senior_key): r.comp.senior_key for r in later}.values()
            multi_groups += [[r for r in later if r.comp.senior_key == key] for key in keys]
        for multi_group in multi_groups:
            multi = [d for d in _multiparent_decompositions(ctx, multi_group) if _valid_decomposition(ctx, d, multi_sets)]
            scored = []
            for parts in _preferred(multi) if multi else []:
                result = _render_multiparent(ctx, parts)
                if result is not None:
                    scored.append((result[0], result[1], parts))
            if scored:
                scored.sort(key=lambda r: (r[0], r[1]))
                _, text, parts = scored[0]
                parents = [p for p in parts if p.role == "parent"]
                return text + _multiparent_name(parents[0], len(parents)), parts[0], ()
        decomps = []
        for root in group:
            decomps.extend(_tree_decompositions(ctx, root))
        decomps = [d for d in decomps if _valid_decomposition(ctx, d, retained_sets)]
        if not decomps:
            continue
        scored = []
        for parts in _preferred(decomps):
            choice = _root_candidates(ctx, parts[0])
            if choice is None:
                continue
            scored.append((choice[0], choice[1].text, parts, choice))
        if not scored:
            continue
        scored.sort(key=lambda s: (s[0], s[1]))
        key, text, parts, choice = scored[0]
        root = parts[0]
        return text + _parent_name(root, True), root, key
    raise UnsupportedStructure("this ring system cannot be named by fusion nomenclature (P-25.5)")


_STANDARD_BONDING = {
    "N": 3, "P": 3, "As": 3, "Sb": 3, "Bi": 3, "B": 3, "Al": 3, "Ga": 3, "In": 3, "Tl": 3,
    "O": 2, "S": 2, "Se": 2, "Te": 2, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4,
    "F": 1, "Cl": 1, "Br": 1, "I": 1,
}


def _kekule(mol):
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Exception:
        pass
    return kekule


def _double_bonds(atom):
    return sum(1 for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 2.0)


def _capacity(atom):
    """Double bonds a ring atom can take in the mancude parent: 1 for a carbon or a divalent-type heteroatom with two
    ring bonds; more for an atom with a nonstandard bonding number (the lambda convention, P-25.6)."""
    symbol = atom.GetSymbol()
    degree = atom.GetDegree()
    if symbol == "C" or symbol in ("Si", "Ge", "Sn", "Pb"):
        return 1
    standard = _STANDARD_BONDING.get(symbol)
    if standard is None:
        return 0
    valence = int(round(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()))
    if valence > standard:
        return valence - degree - atom.GetTotalNumHs() if valence - degree - atom.GetTotalNumHs() > 0 else 0
    return 1 if standard - degree - atom.GetTotalNumHs() == 1 and degree <= 2 else 0


def _valence(atom):
    return int(round(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()))


def _above_standard(atom):
    standard = _STANDARD_BONDING.get(atom.GetSymbol())
    return standard is not None and _valence(atom) > standard


def _bonding_numbers(mol):
    """{atom: n} of ring heteroatoms with a bonding number above the standard one (P-25.6) and {atom: c} of atoms with c
    double bonds in the ring, c >= 2 (the delta convention, P-25.7.2)."""
    lam, delta = {}, {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        valence = int(round(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()))
        if symbol in _STANDARD_BONDING and valence > _STANDARD_BONDING[symbol]:
            lam[atom.GetIdx()] = valence
        if _double_bonds(atom) >= 2:
            delta[atom.GetIdx()] = _double_bonds(atom)
    return lam, delta


def _indicated_hydrogen_atoms(mol):
    """Ring atoms that could carry a double bond in the mancude parent but do not: the positions of indicated hydrogen."""
    kekule = _kekule(mol)
    atoms = []
    for atom in kekule.GetAtoms():
        capacity = _capacity(atom)
        if capacity and _double_bonds(atom) < capacity and (atom.GetTotalNumHs() > 0 or capacity > 1):
            atoms.append(atom.GetIdx())
    for atom in kekule.GetAtoms():
        if not _capacity(atom) and atom.GetTotalNumHs() > 0 and atom.IsInRing() and not _above_standard(atom):
            atoms.append(atom.GetIdx())
    return atoms


def system_numbering_options(ctx, name, root):
    if isinstance(root, _ReplacementRoot):
        return root.numberings
    if len(root.rings) == ctx.n:
        return root.numberings_in_system()
    if name == "cyclopenta[a]phenanthrene":
        images = exception_numbering(ctx.sk, name)
        if images:
            return images
    return fused_numberings(ctx.sk)


def fused_parent_data(mol):
    """(name, numbering options, indicated-hydrogen atoms, lambda map, delta map) of the fused ring system `mol`."""
    ctx = Context(mol)
    if ctx.n < 2:
        raise UnsupportedStructure("a single ring is not a fused ring system")
    name, root = fusion_name(mol)
    kekule = _kekule(mol)
    lam, delta = _bonding_numbers(kekule)
    return name, system_numbering_options(ctx, name, root), _indicated_hydrogen_atoms(kekule), lam, delta


def marked_name(name, numbering, indicated, lam, delta, stem=None):
    """Parent name with its indicated-hydrogen, lambda and delta locants; a lambda or delta mark joins the leading
    heteroatom locant of the fusion name when that locant is cited there (P-25.6)."""
    ih_text = ",".join(f"{t}H" for t in sorted((numbering[a] for a in indicated), key=_locant_key))
    leading = re.match(r"(\d+[a-z]?(?:,\d+[a-z]?)*)-(.+)", name)
    cited = leading.group(1).split(",") if leading else []
    rest = leading.group(2) if leading else name
    loose = []
    for a in sorted(set(lam) | set(delta), key=lambda a: _locant_key(numbering[a])):
        mark = (f"\u03bb{lam[a]}" if a in lam else "") + (f"\u03b4{delta[a]}" if a in delta else "")
        if numbering[a] in cited:
            cited[cited.index(numbering[a])] += mark
        else:
            loose.append(f"{numbering[a]}{mark}")
    body = (",".join(cited) + "-" + rest) if cited else rest
    return "-".join(part for part in (ih_text, ",".join(loose), body) if part)


def name_fused_ring_system(mol):
    """Name with indicated hydrogen, lambda and delta locants of the fused ring system `mol` (one connected ortho- and
    peri-fused system)."""
    name, options, indicated, lam, delta = fused_parent_data(mol)
    if not indicated and not lam and not delta:
        return name

    def key(numbering):
        return (
            sorted(_locant_key(numbering[a]) for a in indicated),
            sorted((-lam[a], _locant_key(numbering[a])) for a in lam),
        )

    best = min(key(n) for n in options)
    chosen = [n for n in options if key(n) == best][0]
    return marked_name(name, chosen, indicated, lam, delta)


def _is_mancude(mol):
    """True when the ring system carries the maximum number of noncumulative double bonds (P-25.7.1.1)."""
    kekule = _kekule(mol)
    capacity = {a.GetIdx(): _capacity(a) for a in kekule.GetAtoms() if _capacity(a)}
    graph = nx.Graph()
    for atom, cap in capacity.items():
        graph.add_nodes_from((atom, k) for k in range(cap))
    for bond in kekule.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in capacity and b in capacity:
            graph.add_edges_from(((a, i), (b, j)) for i in range(capacity[a]) for j in range(capacity[b]))
    matched = len(nx.max_weight_matching(graph, maxcardinality=True))
    minimum = sum(capacity.values()) - 2 * matched
    actual = sum(capacity[a] - min(capacity[a], _double_bonds(kekule.GetAtomWithIdx(a))) for a in capacity)
    return actual == minimum and all(
        a.GetTotalNumHs() <= 1 or a.GetIdx() in _indicated_hydrogen_atoms(kekule) for a in kekule.GetAtoms()
    )


def fused_ring_system_name(mol):
    """Fusion name of a bare mancude fused ring system, or None when `mol` is anything else."""
    if mol.GetNumAtoms() < 5 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    info = mol.GetRingInfo()
    if PREFER_VON_BAEYER.get() and sum(len(ring) >= 5 for ring in info.AtomRings()) < 2:
        return None
    if not 2 <= info.NumRings() <= _MAX_RINGS_IN_SYSTEM or any(info.NumAtomRings(i) == 0 for i in range(mol.GetNumAtoms())):
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
        return None
    try:
        if not _is_mancude(mol):
            return None
        return name_fused_ring_system(mol)
    except UnsupportedStructure:
        return None
