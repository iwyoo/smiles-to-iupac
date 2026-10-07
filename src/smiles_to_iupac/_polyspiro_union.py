"""Polyspiro ring systems of three or more components, at least one polycyclic (P-24.4, P-24.6, P-24.7, P-24.8): the
component names are cited in brackets after 'dispiro', 'trispiro', ..., the locants of each component are primed by
its position in the citation order, and the spiro atoms are cited as locant pairs. Unbranched systems list the
components in order of occurrence from the alphabetically lower terminal; three identical ones take 'dispiroter'. A
central component with terminal ones is cited by P-24.7.1 and P-24.7.2.
"""

from collections import defaultdict
from itertools import permutations, product

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, multiplied_word
from ._fused_numbering import HETERO_RANK
from ._numerals import multiplying_prefix
from ._spiro_union import (
    _A_PREFIX,
    _HALOGENS,
    _primed,
    _bonding_number,
    _capable,
    _component,
    _components,
    _ene_name,
    _hydrogen_choices,
    _lk,
    _name_key,
    _replacement_prefixes,
    _spiro_atom,
)
from ._substituents import format_substituent_prefixes

_MAX_ASSIGNMENTS = 200000


def _structure(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    comps = _components(mol)
    if len(comps) < 3 or max(c["rings"] for c in comps) < 2:
        return None
    membership = defaultdict(list)
    for i, comp in enumerate(comps):
        for a in comp["atoms"]:
            membership[a].append(i)
    spiro = {a: members for a, members in membership.items() if len(members) > 1}
    if any(len(members) > 3 for members in spiro.values()) or sum(len(m) - 1 for m in spiro.values()) != len(comps) - 1:
        return None
    neighbours = defaultdict(set)
    for members in spiro.values():
        for i in members:
            neighbours[i].update(m for m in members if m != i)
    seen, queue = {0}, [0]
    while queue:
        for n in neighbours[queue.pop()]:
            if n not in seen:
                seen.add(n)
                queue.append(n)
    if len(seen) != len(comps):
        return None
    return comps, spiro, neighbours


def _shape(spiro, count):
    atoms_of = {i: [s for s, members in spiro.items() if i in members] for i in range(count)}
    if all(len(m) == 2 for m in spiro.values()) and max(len(a) for a in atoms_of.values()) <= 2:
        return "chain", None
    centres = [i for i in range(count) if all(i in m for m in spiro.values())]
    others_terminal = lambda c: all(len(atoms_of[i]) == 1 for i in range(count) if i != c)
    if len(spiro) == 1:
        return "hub", None
    if len(centres) == 1 and others_terminal(centres[0]):
        return "star", centres[0]
    return None, None


def has_polyspiro_union_shape(mol) -> bool:
    structure = _structure(mol)
    return structure is not None and _shape(structure[1], len(structure[0]))[0] is not None


def _unique(orders):
    result = []
    for order in orders:
        if order not in result:
            result.append(order)
    return result


def _chain_layouts(spiro, names, count):
    atoms_of = {i: [s for s, members in spiro.items() if i in members] for i in range(count)}
    start = next(c for c in range(count) if len(atoms_of[c]) == 1)
    path = [start]
    while len(path) < count:
        atom = next(s for s in atoms_of[path[-1]] if all(m not in path[:-1] for m in spiro[s]) and any(m not in path for m in spiro[s]))
        path.append(next(m for m in spiro[atom] if m not in path))
    keyed = [([_name_key(names[c]) for c in o], o) for o in (path, path[::-1])]
    best = min(k for k, _ in keyed)
    layouts = []
    for order in _unique([o for k, o in keyed if k == best]):
        links = [(next(s for s, m in spiro.items() if {x, y} == set(m)), x, y) for x, y in zip(order, order[1:])]
        layouts.append({"kind": "chain", "flat": order, "groups": [[c] for c in order], "links": links})
    return layouts


def _classes(members, names):
    classes = defaultdict(list)
    for c in members:
        classes[names[c]].append(c)
    return sorted(classes.values(), key=lambda cls: _name_key(names[cls[0]]))


def _star_layouts(spiro, names, centre, count):
    terminals = [c for c in range(count) if c != centre]
    atom_of = {t: next(s for s, m in spiro.items() if t in m) for t in terminals}
    classes = _classes(terminals, names)

    def links_for(flat):
        return [(atom_of[t], centre, t) for t in flat if t != centre]

    layouts = []
    if len(classes) == 1:
        for order in permutations(terminals):
            flat = [centre, *order]
            layouts.append({"kind": "star", "flat": flat, "groups": [[centre], list(order)], "links": links_for(flat), "centre": centre})
        return layouts
    first, rest = classes[0], classes[1:]
    for lead in first:
        mates = [t for t in terminals if t != lead and atom_of[t] == atom_of[lead]]
        if any(t not in first for t in mates):
            raise UnsupportedStructure("different terminal components on one spiro atom are not supported yet")
        remaining_first = [t for t in first if t != lead and t not in mates]
        for mate_order in permutations(mates):
            for first_order in permutations(remaining_first):
                for rest_orders in product(*(permutations(cls) for cls in rest)):
                    flat = [lead, *mate_order, centre, *first_order, *(c for order in rest_orders for c in order)]
                    group_one = [lead, *mate_order, *first_order]
                    groups = [group_one, [centre], *(list(order) for order in rest_orders)]
                    layouts.append({"kind": "star", "flat": flat, "groups": groups, "links": links_for(flat), "centre": centre})
    return layouts


def _hub_layouts(spiro, names, count):
    (atom,) = spiro
    comps = sorted(range(count), key=lambda c: _name_key(names[c]))
    classes = _classes(comps, names)
    if len(classes) == 1:
        return [{"kind": "hub_ter", "flat": comps, "groups": [comps], "links": [], "hub": atom}]
    if len(classes) == 3:
        a, b, c = (cls[0] for cls in classes)
        return [{"kind": "hub_different", "flat": [a, b, c], "groups": [[a], [b], [c]], "links": [(atom, a, b), (atom, a, c)]}]
    pair = next(cls for cls in classes if len(cls) == 2)
    single = next(cls for cls in classes if len(cls) == 1)[0]
    ordered = [pair, [single]] if classes[0] is pair else None
    if ordered is None:
        raise UnsupportedStructure("this arrangement of components on one spiro atom is not supported yet")
    layouts = []
    for order in permutations(pair):
        flat = [*order, single]
        layouts.append({"kind": "hub_pair", "flat": flat, "groups": [list(order), [single]], "links": [(atom, x, single) for x in order]})
    return layouts


def _blocks(layout):
    group_of = {c: gi for gi, group in enumerate(layout["groups"]) for c in group}
    flat_index = {c: i for i, c in enumerate(layout["flat"])}
    blocks = defaultdict(list)
    for s, x, y in layout["links"]:
        early, late = sorted((x, y), key=flat_index.get)
        blocks[max(group_of[early], group_of[late])].append((flat_index[late], flat_index[early], s, early, late))
    return {gi: [entry[2:] for entry in sorted(entries)] for gi, entries in blocks.items()}


def name_polyspiro_union(mol) -> str:
    from ._ylium_ring import _prefixes

    comps, spiro, _ = _structure(mol)
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopes, radicals and stereodescriptors of a spiro union are not supported yet")
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
        raise UnsupportedStructure("double-bond stereo of a spiro union is not supported yet")
    info = {s: _spiro_atom(mol, s) for s in spiro}
    if any(a.GetFormalCharge() and a.GetIdx() not in spiro for a in mol.GetAtoms()):
        raise UnsupportedStructure("a charge away from the spiro atoms is not supported yet")
    cationic = [s for s, (_, charged) in info.items() if charged]
    if len(cationic) > 1:
        raise UnsupportedStructure("several cationic spiro atoms are not supported yet")
    ring_atoms = set().union(*(c["atoms"] for c in comps))
    if not cationic:
        for a in mol.GetAtoms():
            if a.GetIdx() not in ring_atoms and a.GetAtomicNum() != 6 and a.GetAtomicNum() not in _HALOGENS:
                raise UnsupportedStructure("a spiro union with a principal characteristic group is not supported yet")
    graph = adjacency(mol)

    forced = set()
    for s, members in spiro.items():
        if mol.GetAtomWithIdx(s).GetAtomicNum() != 6 and not info[s][0]:
            if not any(comps[m]["von_baeyer"] for m in members):
                raise UnsupportedStructure("a standard-valence heteroatom at the spiro atom is not supported without a bridged component")
            forced.update(members)
    named = [_component(mol, c, i in forced) for i, c in enumerate(comps)]

    inner = defaultdict(list)
    for s, members in spiro.items():
        if mol.GetAtomWithIdx(s).GetAtomicNum() != 6 and info[s][0]:
            replaced = [m for m in members if named[m]["replacement"]]
            if len(replaced) > 1:
                raise UnsupportedStructure("a nonstandard spiro heteroatom between two skeletal replacement components is not supported yet")
            for m in replaced:
                inner[m].append(s)
    names = {}
    for i, n in enumerate(named):
        n["display"] = n["name"]
        if i in inner:
            n["display"] = _replacement_prefixes(mol, inner[i], {a: "1" for a in inner[i]}, with_lambda=False) + n["name"]
        names[i] = n["display"]

    spiro_in = {i: sorted(s for s, members in spiro.items() if i in members) for i in range(len(comps))}
    capable, polycyclic_bonds, monocyclic_bonds, replacement_atoms = set(), set(), set(), set()
    for i, (comp, n) in enumerate(zip(comps, named)):
        if n["fused"]:
            capable |= _capable(mol, comp, set(spiro_in[i]))
            polycyclic_bonds |= comp["bonds"]
        elif n["replacement"]:
            replacement_atoms |= {
                a for a in comp["atoms"] if mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and a not in inner.get(i, ())
            }
        else:
            monocyclic_bonds |= comp["bonds"]
    for s, (lam, charged) in info.items():
        if lam and not charged:
            pi_bonds = lam - mol.GetAtomWithIdx(s).GetDegree()
            if pi_bonds > 1:
                raise UnsupportedStructure("a spiro atom with several multiple bonds is not supported yet")
            if pi_bonds == 1 and any(named[m]["fused"] for m in spiro[s]):
                capable.add(s)
    saturated, choices = _hydrogen_choices(mol, capable, polycyclic_bonds, monocyclic_bonds)

    shape, centre_comp = _shape(spiro, len(comps))
    if shape == "chain":
        layouts = _chain_layouts(spiro, names, len(comps))
    elif shape == "star":
        layouts = _star_layouts(spiro, names, centre_comp, len(comps))
    else:
        layouts = _hub_layouts(spiro, names, len(comps))

    projections = {}
    for i, n in enumerate(named):
        table = defaultdict(list)
        for numbering, ene in n["numberings"]:
            table[tuple(numbering[s] for s in spiro_in[i])].append((numbering, ene))
        projections[i] = table
    ene_comps = [i for i, n in enumerate(named) if n["numberings"][0][1]]

    best = None
    for layout in layouts:
        flat = layout["flat"]
        if any(i != flat[0] for i in ene_comps):
            continue
        prime_of = {c: k for k, c in enumerate(flat)}
        blocks = _blocks(layout)
        centre = layout.get("centre")
        terminal_atoms = [s for s, _, _ in layout["links"]] if centre is not None else []
        indices = list(range(len(comps)))
        shortlist, shortlist_key = [], None
        for combo in product(*(list(projections[i]) for i in indices)):
            base = {}
            for i, proj in zip(indices, combo):
                for s, loc in zip(spiro_in[i], proj):
                    base[(s, i)] = loc
            locant = {k: _primed(loc, prime_of[k[1]]) for k, loc in base.items()}
            central = tuple(_lk(base[(s, centre)]) for s in terminal_atoms) if centre is not None and not any(l for l, _ in info.values()) else ()
            citation = tuple(
                _lk(locant[(s, c)]) for gi in sorted(blocks) for s, early, late in blocks[gi] for c in (early, late)
            )
            key = (central, sorted(_lk(l) for l in locant.values()), citation)
            if shortlist_key is None or key < shortlist_key:
                shortlist, shortlist_key = [combo], key
            elif key == shortlist_key:
                shortlist.append(combo)
        for combo in shortlist:
            tables = [projections[i][proj] for i, proj in zip(indices, combo)]
            size = 1
            for table in tables:
                size *= len(table)
            if size > _MAX_ASSIGNMENTS:
                raise UnsupportedStructure("too many numberings of the spiro components")
            for assignment in product(*tables):
                locant_of = {}
                for i, (numbering, _) in zip(indices, assignment):
                    for a, loc in numbering.items():
                        text = _primed(loc, prime_of[i])
                        if a not in locant_of or _lk(text) < _lk(locant_of[a]):
                            locant_of[a] = text
                hetero = sorted(replacement_atoms, key=lambda a: (HETERO_RANK[mol.GetAtomWithIdx(a).GetSymbol()], _lk(locant_of[a])))
                heteroatoms = (
                    sorted(_lk(locant_of[a]) for a in hetero),
                    [_lk(locant_of[a]) for a in hetero],
                    [-(_bonding_number(mol.GetAtomWithIdx(a)) or 0) for a in sorted(hetero, key=lambda a: _lk(locant_of[a]))],
                )
                grouped = _prefixes(mol, graph, ring_atoms, locant_of)
                locant_set = sorted(_lk(l) for info_ in grouped.values() for l in info_["locants"])
                citation_prefix = tuple(
                    _lk(l) for name in sorted(grouped, key=alpha_sort_key) for l in sorted(grouped[name]["locants"], key=_lk)
                )
                choice = min(
                    (
                        (sorted(_lk(locant_of[a]) for a in picked), sorted(_lk(locant_of[a]) for a in saturated - picked), picked)
                        for picked in choices
                    ),
                    key=lambda c: (c[0], c[1]),
                )
                ene = assignment[flat[0]][1]
                unsaturation = sorted(choice[1] + [_lk(str(x)) for x in ene])
                key = (shortlist_key, heteroatoms, choice[0], unsaturation, locant_set, citation_prefix)
                if best is None or key < best[0]:
                    best = (key, layout, assignment, locant_of, choice[2], grouped, ene, blocks)
    if best is None:
        raise UnsupportedStructure("this unsaturated component cannot be cited first")
    _, layout, assignment, locant_of, picked, grouped, ene, blocks = best
    flat, groups = layout["flat"], layout["groups"]
    prime_of = {c: k for k, c in enumerate(flat)}

    def locant_in(s, c):
        return _primed(assignment[c][0][s], prime_of[c])

    def shown(c):
        text = named[c]["name"]
        if c == flat[0] and ene:
            text = _ene_name(text, ene)
        if c in inner:
            prefix = _replacement_prefixes(mol, inner[c], {a: locant_in(a, c) for a in inner[c]}, with_lambda=False)
            text = prefix + text
        return text

    def group_text(group):
        texts = {shown(c) for c in group}
        if len(texts) > 1:
            raise UnsupportedStructure("identical components with different heteroatom locants are not supported yet")
        text = texts.pop()
        return text if len(group) == 1 else f"{multiplying_prefix(len(group), compound=True)}({text})"

    def pair_text(gi):
        return ":".join(f"{locant_in(s, early)},{locant_in(s, late)}" for s, early, late in blocks[gi])

    indicated = ",".join(f"{locant_of[a]}H" for a in sorted(picked, key=lambda a: _lk(locant_of[a])))
    hydro_atoms = sorted(saturated - picked, key=lambda a: _lk(locant_of[a]))
    hydro = f"{','.join(locant_of[a] for a in hydro_atoms)}-{multiplied_word(len(hydro_atoms), 'hydro')}" if hydro_atoms else ""
    prefixes = format_substituent_prefixes(grouped) if grouped else ""

    lam_front = ",".join(
        f"{lowest}λ{info[s][0]}"
        for lowest, s in sorted(
            ((min((locant_in(s, m) for m in spiro[s]), key=_lk), s) for s in spiro if info[s][0]),
            key=lambda pair: _lk(pair[0]),
        )
    )
    count = multiplying_prefix(len(spiro)) if len(spiro) > 1 else ""
    kind = layout["kind"]
    identical_chain = kind == "chain" and len(comps) == 3 and len({names[c] for c in range(3)}) == 1
    if kind == "hub_ter":
        hub = layout["hub"]
        locants = sorted((locant_in(hub, c) for c in flat), key=_lk)
        lam = info[hub][0]
        last = shown(flat[0])
        body = f"{','.join([f'{locants[0]}λ{lam}' if lam else locants[0], *locants[1:]])}-spiroter[{last}]"
    elif identical_chain:
        if lam_front:
            raise UnsupportedStructure("a nonstandard spiro atom in three identical components is not supported yet")
        pairs = ":".join(pair_text(gi) for gi in sorted(blocks))
        last = shown(flat[0])
        body = f"{pairs}-{count}spiroter[{last}]"
    else:
        parts = []
        for gi, group in enumerate(groups):
            if gi:
                parts.append(pair_text(gi))
            parts.append(group_text(group))
        last = parts[-1]
        if kind == "hub_different":
            parts[2] = f"({parts[2]})"
        body = f"{count}spiro[{'-'.join(parts)}]"
        if lam_front and kind != "hub_ter":
            body = f"{lam_front}-{body}"
    replacement_prefix = _replacement_prefixes(mol, replacement_atoms, locant_of) if replacement_atoms else ""
    if replacement_prefix:
        body = replacement_prefix + ("-" if body[0].isdigit() else "") + body
    if cationic:
        if last.endswith("e") and not last.endswith(")"):
            body = body[: -len(last) - 1] + last[:-1] + "]"
        lowest = min((locant_in(cationic[0], m) for m in spiro[cationic[0]]), key=_lk)
        body += f"-{lowest}-ylium"

    out = prefixes
    for token in (hydro, indicated, body):
        if not token:
            continue
        joined = out and (token[0].isdigit() or token[0] in "[(" or out.endswith("H"))
        out += ("-" if joined else "") + token
    return out
