"""Metallocenes as parent hydrides (P-69.2.7): the retained -ocene names of
Fe/Ru/Os/Ni/Cr/Co/V sandwiches, substituted by prefixes and suffixes with the
ring atoms numbered 1-5 and 1'-5' ('1,1'-dimethylferrocene', 'ferrocen-1-ol',
'ferrocene-1,1'-dicarboxylic acid'), or cited as a substituent group
('1'-methylferrocen-1-yl') on a parent named by the other modules.
"""

import itertools

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._numerals import multiplying_prefix, numerical_term
from ._phosphanyl_group import PREFIX_PROP
from ._pin import mark
from ._prefix_groups import PrefixNamer, enclose
from ._substituents import format_substituent_prefixes

OCENES = {
    "Fe": "ferrocene", "Ru": "ruthenocene", "Os": "osmocene", "Ni": "nickelocene",
    "Cr": "chromocene", "Co": "cobaltocene", "V": "vanadocene",
}
_PRIME = "′"
_SUFFIX_ORDER = ["acid", "sulfonic", "amide", "nitrile", "al", "ol", "thiol", "amine"]
_SUFFIX_WORD = {
    "acid": "carboxylic acid", "sulfonic": "sulfonic acid", "amide": "carboxamide", "nitrile": "carbonitrile",
    "al": "carbaldehyde", "ol": "ol", "thiol": "thiol", "amine": "amine",
}
_PREFIX_WORD = {
    "acid": "carboxy", "sulfonic": "sulfo", "amide": "carbamoyl", "nitrile": "cyano", "al": "formyl",
    "ol": "hydroxy", "thiol": "sulfanyl", "amine": "amino",
}


def _lk(locant):
    return locant


def _fmt(locant):
    return f"{locant[0]}{_PRIME * locant[1]}"


def _is_cp(mol, ring):
    ring = set(ring)
    if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return True
    doubles = sum(
        1 for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring and b.GetBondTypeAsDouble() == 2.0
    )
    return doubles >= 2


def _cp_rings(mol, graph):
    rings = []
    info = mol.GetRingInfo()
    for ring in info.AtomRings():
        if len(ring) != 5 or any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in ring):
            continue
        if any(info.NumAtomRings(a) > 1 for a in ring) or not _is_cp(mol, ring):
            continue
        rings.append(ring)
    return rings


def _balanced_cyclopentadienides(mol, rings):
    """Two cyclopentadienyl rings (anions or radicals) and a metal ion that makes the whole neutral are the ocene
    (P-69.2.7); any other charge or a neutral cyclopentadiene is an ionic or adduct name (P-14.8.2)."""
    ring_valence = lambda ring: sum(mol.GetAtomWithIdx(a).GetNumRadicalElectrons() - mol.GetAtomWithIdx(a).GetFormalCharge() for a in ring)
    return sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 0 and all(ring_valence(ring) == 1 for ring in rings)


def find_ocene(mol):
    metals = [a for a in mol.GetAtoms() if a.GetSymbol() in OCENES]
    if len(metals) != 1:
        return None
    metal = metals[0]
    graph = adjacency(mol)
    rings = _cp_rings(mol, graph)
    bonded = [r for r in rings if all(metal.GetIdx() in graph[a] for a in r)]
    if len(bonded) == 2 and metal.GetDegree() == 10:
        chosen = bonded
    elif metal.GetDegree() == 0 and len(rings) == 2 and _balanced_cyclopentadienides(mol, rings):
        chosen = rings
    else:
        return None
    if len(Chem.GetMolFrags(mol)) > 1 and metal.GetDegree() != 0:
        return None
    return metal.GetIdx(), [list(r) for r in chosen], graph


def _paired_ocenes(mol):
    metals = [a for a in mol.GetAtoms() if a.GetSymbol() in OCENES]
    if len(metals) < 2:
        return None
    graph = adjacency(mol)
    rings = _cp_rings(mol, graph)
    units = []
    for metal in metals:
        bonded = [list(r) for r in rings if all(metal.GetIdx() in graph[a] for a in r)]
        if len(bonded) != 2 or metal.GetDegree() != 10:
            return None
        units.append((metal.GetIdx(), bonded))
    return units, graph


def _benzo_ocene(mol):
    """A cyclopentadienyl ring plus an indenyl ring on one ocene metal."""
    metals = [a for a in mol.GetAtoms() if a.GetSymbol() in OCENES]
    if len(metals) != 1 or mol.GetNumAtoms() != 15 or metals[0].GetDegree() != 10:
        return None
    metal = metals[0]
    info = mol.GetRingInfo()
    five = [r for r in info.AtomRings() if len(r) == 5 and all(metal.GetIdx() not in (a,) for a in r)]
    bonded = [r for r in five if all(mol.GetBondBetweenAtoms(metal.GetIdx(), a) for a in r)]
    carbons = [a for a in mol.GetAtoms() if a.GetIdx() != metal.GetIdx()]
    if len(bonded) != 2 or any(a.GetAtomicNum() != 6 for a in carbons):
        return None
    fused = [r for r in bonded if any(info.NumAtomRings(a) > 1 for a in r)]
    benzene = [r for r in info.AtomRings() if len(r) == 6]
    if len(fused) != 1 or len(benzene) != 1 or set(benzene[0]) & set(range(15)) != set(benzene[0]):
        return None
    return OCENES[metal.GetSymbol()]


def has_ocene_shape(mol) -> bool:
    return find_ocene(mol) is not None or _paired_ocenes(mol) is not None or _benzo_ocene(mol) is not None


def _cycle_order(graph, ring, start, direction):
    ring_set = set(ring)
    order, prev = [start], None
    while len(order) < 5:
        nbrs = [v for v in graph[order[-1]] if v in ring_set and v != prev]
        nxt = nbrs[direction] if len(nbrs) > 1 else nbrs[0]
        prev = order[-1]
        order.append(nxt)
    return order


def _numberings(graph, rings):
    for first, second in itertools.permutations(range(2)):
        for s1 in rings[first]:
            for d1 in (0, 1):
                for s2 in rings[second]:
                    for d2 in (0, 1):
                        positions = {}
                        for k, a in enumerate(_cycle_order(graph, rings[first], s1, d1), start=1):
                            positions[a] = (k, 0)
                        for k, a in enumerate(_cycle_order(graph, rings[second], s2, d2), start=1):
                            positions[a] = (k, 1)
                        yield positions


def _suffix_kind(mol, graph, ring_atom, root):
    atom = mol.GetAtomWithIdx(root)
    z = atom.GetAtomicNum()
    others = [x for x in graph[root] if x != ring_atom]
    if z == 8 and atom.GetTotalNumHs() == 1 and not others:
        return "ol"
    if z == 16 and atom.GetTotalNumHs() == 1 and not others:
        return "thiol"
    if z == 7 and atom.GetTotalNumHs() == 2 and not others and not atom.GetFormalCharge():
        return "amine"
    if z == 6 and not atom.GetIsAromatic() and all(mol.GetAtomWithIdx(x).GetAtomicNum() != 6 for x in others):
        oxo = [x for x in others if mol.GetAtomWithIdx(x).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(root, x).GetBondTypeAsDouble() == 2.0]
        triple_n = [x for x in others if mol.GetAtomWithIdx(x).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(root, x).GetBondTypeAsDouble() == 3.0]
        if triple_n and len(others) == 1:
            return "nitrile"
        if len(oxo) == 1:
            rest = [x for x in others if x != oxo[0]]
            if not rest and atom.GetTotalNumHs() == 1:
                return "al"
            if len(rest) == 1:
                r = mol.GetAtomWithIdx(rest[0])
                if r.GetAtomicNum() == 8 and r.GetTotalNumHs() == 1 and r.GetDegree() == 1:
                    return "acid"
                if r.GetAtomicNum() == 7 and r.GetTotalNumHs() == 2 and r.GetDegree() == 1:
                    return "amide"
    return None


def _components(graph, unit, roots):
    """{root: atoms of the substituent component reached from root}."""
    out = {}
    for ring_atom, root in roots:
        seen, stack = {root}, [root]
        while stack:
            for v in graph[stack.pop()]:
                if v not in unit and v not in seen:
                    seen.add(v)
                    stack.append(v)
        out[root] = seen
    return out


def _senior_elsewhere(mol, graph, unit, ring_roots, kind):
    """True if a characteristic group more senior than `kind` exists outside
    the ring-attached suffix groups (so the ocene is not the parent)."""
    rank = _SUFFIX_ORDER.index(kind)
    skip = set()
    for _, root in ring_roots:
        skip.add(root)
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        if i in unit or i in skip:
            continue
        z = atom.GetAtomicNum()
        if z == 6 and not atom.GetIsAromatic():
            oxo = [x for x in graph[i] if mol.GetAtomWithIdx(x).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(i, x).GetBondTypeAsDouble() == 2.0]
            hetero = [x for x in graph[i] if mol.GetAtomWithIdx(x).GetAtomicNum() in (7, 8) and x not in oxo]
            if oxo:
                has_oh = any(mol.GetAtomWithIdx(x).GetAtomicNum() == 8 and mol.GetAtomWithIdx(x).GetTotalNumHs() == 1 for x in hetero)
                if has_oh and _SUFFIX_ORDER.index("acid") < rank:
                    return True
                if any(mol.GetAtomWithIdx(x).GetAtomicNum() == 7 for x in hetero) and _SUFFIX_ORDER.index("amide") < rank:
                    return True
                if mol.GetAtomWithIdx(i).GetTotalNumHs() == 1 and not hetero and _SUFFIX_ORDER.index("al") < rank:
                    return True
        if z == 6 and any(mol.GetAtomWithIdx(x).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(i, x).GetBondTypeAsDouble() == 3.0 for x in graph[i]):
            if _SUFFIX_ORDER.index("nitrile") < rank:
                return True
    return False


def _suffix_text(stem, kind, locants):
    word = _SUFFIX_WORD[kind]
    locs = ",".join(_fmt(l) for l in sorted(locants))
    if len(locants) == 1:
        if kind in ("ol", "amine", "thiol"):
            return f"{stem[:-1]}-{locs}-{word}"
        return f"{stem}-{locs}-{word}"
    mult = multiplying_prefix(len(locants))
    return f"{stem}-{locs}-{mult}{word}"


_ORDER = ["acid", "sulfonic", "ester", "amide", "nitrile", "al", "one", "ol", "thiol", "amine"]


def _group_kinds(mol, graph, atoms):
    """One characteristic-group class (suffix candidate) per group found in the atoms."""
    kinds = []
    for i in atoms:
        atom = mol.GetAtomWithIdx(i)
        z = atom.GetAtomicNum()
        kind = None
        if z == 6 and not atom.GetIsAromatic():
            nbrs = [x for x in graph[i]]
            oxo = [x for x in nbrs if mol.GetAtomWithIdx(x).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(i, x).GetBondTypeAsDouble() == 2.0]
            if any(mol.GetAtomWithIdx(x).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(i, x).GetBondTypeAsDouble() == 3.0 for x in nbrs):
                kind = "nitrile"
            elif oxo:
                hetero = [x for x in nbrs if x not in oxo and mol.GetAtomWithIdx(x).GetAtomicNum() in (7, 8, 16)]
                if hetero:
                    h = mol.GetAtomWithIdx(hetero[0])
                    if h.GetAtomicNum() == 8:
                        kind = "acid" if h.GetTotalNumHs() == 1 else "ester"
                    elif h.GetAtomicNum() == 7:
                        kind = "amide"
                elif atom.GetTotalNumHs() == 1:
                    kind = "al"
                else:
                    kind = "one"
        elif z == 8 and atom.GetTotalNumHs() == 1 and graph[i]:
            nb = mol.GetAtomWithIdx(graph[i][0])
            if nb.GetAtomicNum() == 6 and not any(mol.GetAtomWithIdx(x).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(nb.GetIdx(), x).GetBondTypeAsDouble() == 2.0 for x in graph[nb.GetIdx()]):
                kind = "ol"
        elif z == 16 and atom.GetTotalNumHs() == 1:
            kind = "thiol"
        elif z == 7 and not atom.GetFormalCharge() and not atom.GetIsAromatic():
            nbrs = [mol.GetAtomWithIdx(x) for x in graph[i]]
            if nbrs and all(n.GetAtomicNum() == 6 and not any(
                mol.GetAtomWithIdx(y).GetAtomicNum() in (7, 8) and mol.GetBondBetweenAtoms(n.GetIdx(), y).GetBondTypeAsDouble() >= 2.0
                for y in graph[n.GetIdx()]
            ) for n in nbrs):
                kind = "amine"
        if kind:
            kinds.append(kind)
    return kinds


def _component_kind(mol, graph, atoms):
    return min(_group_kinds(mol, graph, atoms), key=_ORDER.index, default=None)


def _chain_path(graph, comp, start, goal):
    parents, stack = {start: None}, [start]
    while stack:
        u = stack.pop()
        for v in graph[u]:
            if v in comp and v not in parents:
                parents[v] = u
                stack.append(v)
    if goal not in parents:
        return None
    path, node = [], goal
    while node is not None:
        path.append(node)
        node = parents[node]
    return path[::-1]


def _ring_locants(graph, ring, link_atom, direction, primes):
    return {a: (k, primes) for k, a in enumerate(_cycle_order(graph, ring, link_atom, direction), start=1)}


def _name_phane(mol, units, graph, parent):
    """Ocenes joined 1,1'- in one ring through carbon nodes (P-52.2.5, P-26):
    '1,3(1,1')-diferrocenacyclotetraphane'; substituent locants on amplificant
    atoms are composite ('1^2')."""
    unit_atoms = [{m} | set(r[0]) | set(r[1]) for m, r in units]
    every = set().union(*unit_atoms)
    owner = {a: k for k, atoms in enumerate(unit_atoms) for a in atoms}
    ring_of = {a: ri for _, rs in units for ri, r in enumerate(rs) for a in r}
    roots = [(a, x) for m, rs in units for r in rs for a in r for x in graph[a] if x not in every and x != m]
    comps = _components(graph, every, roots)
    by_comp = {}
    for a, x in roots:
        by_comp.setdefault(frozenset(comps[x]), []).append((a, x))
    links, ring_subs = [], []
    for comp, ends in by_comp.items():
        if len(ends) == 2 and owner[ends[0][0]] != owner[ends[1][0]]:
            links.append((ends, comp))
        elif len(ends) == 1:
            ring_subs.append(ends[0])
        else:
            return None
    count = len(units)
    if len(links) != count:
        return None
    namer = PrefixNamer(mol, graph)
    link_at = {}
    chains = []
    chain_subs = {}
    for li, (ends, comp) in enumerate(links):
        for side, (a, x) in enumerate(ends):
            if (owner[a], ring_of[a]) in link_at:
                return None
            link_at[(owner[a], ring_of[a])] = (li, side)
        path = _chain_path(graph, comp, ends[0][1], ends[1][1])
        if path is None or any(
            mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).GetIsAromatic() for i in path
        ) or any(b.GetBondTypeAsDouble() != 1.0 for i in comp for b in mol.GetAtomWithIdx(i).GetBonds()):
            return None
        chains.append(path)
        for c in path:
            for n in graph[c]:
                if n in comp and n not in path and n not in chain_subs.get(c, []):
                    if _component_kind(mol, graph, _group_atoms(graph, n, set(path) | every)):
                        return None
                    chain_subs.setdefault(c, []).append(n)
    if len(link_at) != 2 * count:
        return None
    for a, x in ring_subs:
        if _suffix_kind(mol, graph, a, x) or _component_kind(mol, graph, comps[x]):
            return None

    seq, unit, ring = [], 0, 0
    while True:
        li, side = link_at[(unit, ring)]
        path = chains[li] if side == 0 else chains[li][::-1]
        seq.append((unit, ring, path))
        (a, _), = [links[li][0][1 - side]]
        unit, ring = owner[a], 1 - ring_of[a]
        if unit == 0:
            break
    if len(seq) != count or ring != 0:
        return None

    forward = seq
    backward = [(seq[i][0], 1 - seq[i][1], seq[i - 1][2][::-1]) for i in range(count - 1, -1, -1)]
    best = None
    for stations in (forward, backward):
        for start in range(count):
            order = stations[start:] + stations[:start]
            subs, positions, pos = [], [], 1
            node_of = {}
            for u, exit_ring, path in order:
                positions.append(pos)
                node_of[u] = pos
                for off, c in enumerate(path, start=1):
                    for n in chain_subs.get(c, []):
                        subs.append((pos + off, namer.name(n, c)))
                pos += 1 + len(path)
            unit_choice = {}
            for u, exit_ring, path in order:
                rs = units[u][1]
                link_atoms = {r: next(a for a, _ in [e for ends, _ in links for e in ends] if owner[a] == u and ring_of[a] == r) for r in (0, 1)}
                local = [(a, x) for a, x in ring_subs if owner[a] == u]
                options = []
                for unprimed in (0, 1):
                    for d0 in (0, 1):
                        for d1 in (0, 1):
                            loc = {}
                            loc.update(_ring_locants(graph, rs[unprimed], link_atoms[unprimed], d0, 0))
                            loc.update(_ring_locants(graph, rs[1 - unprimed], link_atoms[1 - unprimed], d1, 1))
                            named = [(loc[a], namer.name(x, a)) for a, x in local]
                            options.append((sorted((l, n[0]) for l, n in named), named))
                unit_choice[u] = min(options, key=lambda o: o[0])
            ring_entries = [(node_of[u], l, n) for u, (_, named) in unit_choice.items() for l, n in named]
            key = (
                positions,
                sorted([p[0] for p in subs] + [e[0] for e in ring_entries]),
                sorted([(p[0], 0, 0) for p in subs] + [(e[0], e[1][0], e[1][1]) for e in ring_entries]),
            )
            if best is None or key < best[0]:
                best = (key, positions, subs, ring_entries)
    _, positions, subs, ring_entries = best
    grouped = {}
    for loc, (name, compound) in subs:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(str(loc))
    for node, (k, primes), (name, compound) in ring_entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(f"{node}^{k}{_PRIME * primes}")
    for entry in grouped.values():
        entry["locants"].sort(key=lambda t: (int(t.split("^")[0]), t))
    total = sum(1 + len(p) for _, _, p in forward)
    prefixes = format_substituent_prefixes(grouped) if grouped else ""
    mult = multiplying_prefix(count)
    stem = parent[:-1] + "a"
    return f"{prefixes}{'-' if prefixes else ''}{','.join(map(str, positions))}(1,1{_PRIME})-{mult}{stem}cyclo{numerical_term(total)}phane"


def _name_paired_ocenes(mol, units, graph):
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    parents = {OCENES[mol.GetAtomWithIdx(m).GetSymbol()] for m, _ in units}
    if len(parents) == 1:
        phane = _name_phane(mol, units, graph, parents.pop())
        if phane:
            return phane
    from ._ocene_linked import name_assembly

    return name_assembly(mol, units, graph)


def name_ocene(mol) -> str:
    benzo = _benzo_ocene(mol)
    if benzo:
        return "benzo" + benzo
    paired = _paired_ocenes(mol)
    if paired is not None:
        return _name_paired_ocenes(mol, *paired)
    metal_idx, rings, graph = find_ocene(mol)
    metal = mol.GetAtomWithIdx(metal_idx)
    parent = OCENES[metal.GetSymbol()]
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    unit = {metal_idx} | set(rings[0]) | set(rings[1])
    roots = [(a, x) for r in rings for a in r for x in graph[a] if x not in unit and x != metal_idx]
    if not roots:
        return parent
    comps = _components(graph, unit, roots)
    bridge = _bridge_roots(roots, comps)
    if bridge:
        return _as_bridge(mol, graph, rings, unit, roots, comps, bridge, parent)
    suffix = {(a, x): _suffix_kind(mol, graph, a, x) for a, x in roots}
    comp_kinds = {(a, x): _group_kinds(mol, graph, comps[x]) for a, x in roots}
    ring_kinds = [k for k in suffix.values() if k]
    ring_best = min(ring_kinds, key=_SUFFIX_ORDER.index, default=None)
    ring_rank = _ORDER.index(ring_best) if ring_best else len(_ORDER)
    chain_roots = [(a, x) for a, x in roots if not suffix[(a, x)] and comp_kinds[(a, x)]]
    chain_best = min((k for r in chain_roots for k in comp_kinds[r]), key=_ORDER.index, default=None)
    chain_rank = _ORDER.index(chain_best) if chain_best else len(_ORDER)
    chain_count = sum(comp_kinds[r].count(chain_best) for r in chain_roots)
    namer = PrefixNamer(mol, graph)

    if chain_best and (chain_rank < ring_rank or (chain_rank == ring_rank and chain_count > ring_kinds.count(ring_best))):
        holders = [r for r in chain_roots if chain_best in comp_kinds[r]]
        try:
            if len(holders) == 2 and _identical_arms(mol, comps, holders):
                others = {r: _root_name(namer, suffix, r) for r in roots if r not in holders}
                return _as_multiplied_parent(mol, graph, rings, unit, comps, holders, others, parent)
            for pick in sorted(holders, key=lambda r: (-comp_kinds[r].count(chain_best), -len(comps[r[1]]))):
                named = {r: _root_name(namer, suffix, r) for r in roots if r != pick}
                try:
                    return _as_substituent_group(mol, graph, rings, unit, named, pick, parent)
                except UnsupportedStructure:
                    continue
        except UnsupportedStructure:
            pass
        named = {r: _root_name(namer, suffix, r) for r in roots}
        return mark(
            _as_parent(mol, graph, rings, roots, suffix, None, named, parent),
            "the senior group of a substituent is cited as a prefix because the substituent could not be named as the parent (P-41)",
        )
    named = {r: _root_name(namer, suffix, r) for r in roots if not (ring_best and suffix[r] == ring_best)}
    return _as_parent(mol, graph, rings, roots, suffix, ring_best, named, parent)


def _root_name(namer, suffix, root):
    kind = suffix[root]
    return (_PREFIX_WORD[kind], False) if kind else namer.name(root[1], root[0])


def _bridge_roots(roots, comps):
    """The two roots joined by one substituent component, a chain or ring bridging the two cyclopentadienyl rings."""
    seen = {}
    for a, x in roots:
        seen.setdefault(frozenset(comps[x]), []).append((a, x))
    bridges = [ends for ends in seen.values() if len(ends) > 1]
    if not bridges:
        return None
    if len(bridges) > 1 or len(bridges[0]) != 2 or bridges[0][0][0] == bridges[0][1][0]:
        raise UnsupportedStructure("several bridges between the rings of one metallocene are not supported here")
    return bridges[0]


def _identical_arms(mol, comps, holders):
    keys = {_component_smiles(mol, comps[x] | {a}, a) for a, x in holders}
    return len(keys) == 1


def _grouped(positions, named):
    grouped: dict = {}
    for (a, x), (name, compound) in named.items():
        entry = grouped.setdefault(name, {"locants": [], "compound": compound})
        entry["locants"].append(_fmt(positions[a]))
    for entry in grouped.values():
        entry["locants"].sort()
    return grouped


def _as_parent(mol, graph, rings, roots, kinds, principal, named, parent):
    best = None
    for positions in _numberings(graph, rings):
        suffix = sorted(positions[a] for a, x in roots if principal and kinds[(a, x)] == principal)
        subs = sorted(positions[a] for (a, x) in named)
        grouped = _grouped(positions, named)
        key = (suffix, subs, [grouped[k]["locants"] for k in sorted(grouped)])
        if best is None or key < best[0]:
            best = (key, suffix, grouped)
    _, suffix, grouped = best
    prefixes = format_substituent_prefixes(grouped)
    stem = _suffix_text(parent, principal, suffix) if principal else parent
    return f"{prefixes}{stem}"


def _stand_in(rw, atom_idx, group):
    """Attach a stand-in for `group` to the atom: a placeholder prefix, or a phenyl ring where a placeholder would read
    as a halogen on a heteroatom or an acyl carbon."""
    atom = rw.GetAtomWithIdx(atom_idx)
    needs_carbon = atom.GetAtomicNum() != 6 or any(
        b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() != 6 for b in atom.GetBonds()
    )
    if not needs_carbon:
        marker = rw.AddAtom(Chem.Atom(53))
        rw.AddBond(atom_idx, marker, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(marker).SetProp(PREFIX_PROP, enclose(group))
        return False
    ring = [rw.AddAtom(Chem.Atom(6)) for _ in range(6)]
    for k, idx in enumerate(ring):
        rw.GetAtomWithIdx(idx).SetIsAromatic(True)
        rw.AddBond(idx, ring[(k + 1) % 6], Chem.BondType.AROMATIC)
    rw.AddBond(atom_idx, ring[0], Chem.BondType.SINGLE)
    return True


def swap_stand_ins(name, stands, once=False):
    """Put the groups for the phenyl rings that stood in for them (see `_stand_in`)."""
    import re

    if not stands:
        return name
    found = list(re.finditer(r"(di|tri|tetra|penta|hexa)?phenyl", name))
    if len(set(stands)) != 1 or len(found) != 1 or "phenylene" in name or bool(found[0].group(1)) != (len(stands) > 1 and not once):
        raise UnsupportedStructure("the metallocene groups could not be cited as prefixes here")
    word = "bis" if len(stands) > 1 and not once else ""
    return f"{name[:found[0].start()]}{word}{enclose(stands[0])}{name[found[0].end():]}"


def _name_contracted(rw, drop, marks, once=False):
    """Name of the molecule left after `drop` is removed, each (atom, group) of `marks` standing for the group; with
    `once` a divalent group is cited a single time at the locants of its two marks."""
    import re

    from .core import _name_mol

    stands = [group for atom_idx, group in marks if _stand_in(rw, atom_idx, group)]
    for idx in sorted(drop, reverse=True):
        rw.RemoveAtom(idx)
    contracted = rw.GetMol()
    Chem.SanitizeMol(contracted)
    name = swap_stand_ins(_name_mol(contracted), stands, once)
    groups = {group for _, group in marks}
    if once:
        group = next(iter(groups))
        match = re.search(r"(\d+(?:,\d+)*)-(?:bis|di)?[(\[]*" + re.escape(group) + r"[)\]]*", name)
        if match:
            name = f"{name[:match.start()]}{match.group(1)}-({group}){name[match.end():]}"
    if "iodo" in name or "iodide" in name or any(group not in name for group in groups):
        raise UnsupportedStructure("the metallocene groups could not be cited as prefixes here")
    return name


def _contract(mol, graph, unit, named, main_root, group_name):
    remove = set(unit)
    for (a, x) in named:
        remove |= _group_atoms(graph, x, unit)
    return _name_contracted(Chem.RWMol(mol), remove, [(main_root, group_name)])


def _as_substituent_group(mol, graph, rings, unit, named, main, parent):
    """The metallocene as a substituent group of the rest of the molecule."""
    main_atom, main_root = main
    best = None
    for positions in _numberings(graph, rings):
        if positions[main_atom] != (1, 0):
            continue
        subs = sorted(positions[a] for (a, x) in named)
        grouped = _grouped(positions, named)
        key = (subs, [grouped[k]["locants"] for k in sorted(grouped)])
        if best is None or key < best[0]:
            best = (key, grouped)
    _, grouped = best
    stem = parent[:-1] if parent.endswith("e") else parent
    group_name = f"{format_substituent_prefixes(grouped)}{stem}-1-yl"
    return _contract(mol, graph, unit, named, main_root, group_name)


def _divalent_group(graph, rings, ends, named, parent):
    """The ocene-1,1'-diyl group with `ends` as its free valences (P-29.2) and the other substituents as prefixes."""
    (a1, _), (a2, _) = ends
    best = None
    for positions in _numberings(graph, rings):
        if positions[a1][0] != 1 or positions[a2][0] != 1 or positions[a1][1] == positions[a2][1]:
            continue
        grouped = _grouped(positions, named)
        key = (sorted(positions[a] for a, _ in named), [grouped[n]["locants"] for n in sorted(grouped)])
        if best is None or key < best[0]:
            best = (key, grouped)
    prefixes = format_substituent_prefixes(best[1]) if best[1] else ""
    return f"{prefixes}{parent}-1,1{_PRIME}-diyl"


def _as_multiplied_parent(mol, graph, rings, unit, comps, holders, others, parent):
    """Two identical characteristic-group components joined by the metallocene
    ('1,1'-(ferrocene-1,1'-diyl)di(ethan-1-one)', P-15.3)."""
    (a1, x1), (a2, x2) = holders
    acyl = any(_is_acyl_oxo(mol, x1, b) for b in mol.GetAtomWithIdx(x1).GetBonds())
    rw = Chem.RWMol(mol)
    remove = set(unit) | comps[x2]
    for r in others:
        remove |= comps[r[1]]
    if acyl:
        ring = [rw.AddAtom(Chem.Atom(6)) for _ in range(6)]
        for k, idx in enumerate(ring):
            rw.GetAtomWithIdx(idx).SetIsAromatic(True)
            rw.AddBond(idx, ring[(k + 1) % 6], Chem.BondType.AROMATIC)
        rw.AddBond(x1, ring[0], Chem.BondType.SINGLE)
        token = "phenyl"
    else:
        token = "QQQ"
        marker = rw.AddAtom(Chem.Atom(53))
        rw.AddBond(x1, marker, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(marker).SetProp(PREFIX_PROP, token)
    for idx in sorted(remove, reverse=True):
        rw.RemoveAtom(idx)
    single = rw.GetMol()
    Chem.SanitizeMol(single)
    from .core import _name_mol
    import re

    name = _name_mol(single)
    match = re.search(r"(\d+)-\(?" + token + r"\)?|\(?" + token + r"\)?", name)
    if not match or name.count(token) != 1:
        raise UnsupportedStructure("the multiplied parent could not be isolated")
    locant = match.group(1)
    stripped = (name[: match.start()] + name[match.end():]).strip("-")
    if locant is None and stripped.startswith("acet"):
        locant = "2"
    inner = f"({stripped})" if acyl or any(ch.isdigit() for ch in stripped) else stripped
    divalent = _divalent_group(graph, rings, holders, others, parent)
    locs = f"{locant},{locant}{_PRIME}" if locant else ""
    head = f"{locs}-" if locs else ""
    return f"{head}({divalent})di{inner}"


def _as_bridge(mol, graph, rings, unit, roots, comps, ends, parent):
    """A chain or ring joining both cyclopentadienyl rings, cited as a divalent metallocene group
    ('3,5-(ferrocene-1,1'-diyl)pentanoic acid')."""
    namer = PrefixNamer(mol, graph)
    suffix = {r: _suffix_kind(mol, graph, *r) for r in roots}
    others = {}
    for r in roots:
        if r in ends:
            continue
        if not suffix[r] and _component_kind(mol, graph, comps[r[1]]):
            raise UnsupportedStructure("a bridged metallocene beside a further group of the same rank is not supported here")
        others[r] = _root_name(namer, suffix, r)
    group = _divalent_group(graph, rings, ends, others, parent)
    remove = set(unit)
    for a, x in roots:
        if (a, x) not in ends:
            remove |= comps[x]
    name = _name_contracted(Chem.RWMol(mol), remove, [(x, group) for _, x in ends], once=True)
    if f"-({group})" not in name:
        raise UnsupportedStructure("the name of the bridge does not cite the metallocene group once")
    if not any(_group_kinds(mol, graph, comps[x]) for _, x in ends) and not any(mol.GetAtomWithIdx(i).IsInRing() for _, x in ends for i in comps[x]):
        return mark(name, "a phane name needs two or more ring systems (P-52.2.5.1), so no PIN is defined for a single bridged metallocene")
    return name


def _is_acyl_oxo(mol, atom_idx, bond):
    other = bond.GetOtherAtom(mol.GetAtomWithIdx(atom_idx))
    return bond.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() == 8 and other.GetDegree() == 1


def _component_smiles(mol, atoms, attach):
    rw = Chem.RWMol(mol)
    rw.GetAtomWithIdx(attach).SetAtomMapNum(1)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        rw.RemoveAtom(idx)
    sub = rw.GetMol()
    sub.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(sub)


def _group_atoms(graph, root, unit):
    seen, stack = {root}, [root]
    while stack:
        for v in graph[stack.pop()]:
            if v not in unit and v not in seen:
                seen.add(v)
                stack.append(v)
    return seen
