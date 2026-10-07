"""Spiro unions of rings with at least one polycyclic component (P-24.3 to P-24.8): the component names are cited in
brackets after 'spiro', 'dispiro', ..., each component's locants are primed by its place in the citation order, and the
spiro atoms are cited as locant pairs. Chains, a central component with terminals, atoms shared by three components and
monocyclic units are laid out per the clauses; the same search numbers the system for suffix groups (P-31.1.4).
"""

import re
from collections import defaultdict
from itertools import permutations, product
from types import SimpleNamespace

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, multiplied_word
from ._fused_numbering import HETERO_RANK
from ._numerals import alkane_name, multiplying_prefix, numerical_term
from ._polyspiro import _arc_choice_options, _build_sequence, _chain_direction_candidates
from ._spiro_union import (
    _HALOGENS,
    SpiroLocant,
    _bonding_number,
    _capable,
    _component,
    _components,
    _hydrogen_choices,
    _lk,
    _name_key,
    _primed,
    _replacement_prefixes,
    _spiro_atom,
)
from ._substituents import format_substituent_prefixes

_MAX_ASSIGNMENTS = 200000


def _raw_structure(mol, atoms):
    if atoms is None and len(Chem.GetMolFrags(mol)) != 1:
        return None
    comps = [c for c in _components(mol) if atoms is None or c["atoms"] <= atoms]
    if len(comps) < 2 or max(c["rings"] for c in comps) < 2:
        return None
    membership = defaultdict(list)
    for i, comp in enumerate(comps):
        for a in comp["atoms"]:
            membership[a].append(i)
    spiro = {a: members for a, members in membership.items() if len(members) > 1}
    if any(len(members) > 3 for members in spiro.values()) or sum(len(m) - 1 for m in spiro.values()) != len(comps) - 1:
        return None
    return comps, spiro


def _connected(spiro, count):
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
    return len(seen) == count


def _merge_units(comps, spiro):
    parent = list(range(len(comps)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    internal = {
        s for s, members in spiro.items() if len(members) == 2 and all(comps[m]["rings"] == 1 for m in members)
    }
    for s in internal:
        i, j = spiro[s]
        parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(len(comps)):
        groups[find(i)].append(i)
    if all(len(g) == 1 for g in groups.values()):
        return None
    merged, index_of = [], {}
    for g in groups.values():
        if len(g) == 1:
            merged.append(comps[g[0]])
        else:
            merged.append(
                {
                    "atoms": set().union(*(comps[i]["atoms"] for i in g)),
                    "bonds": set().union(*(comps[i]["bonds"] for i in g)),
                    "rings": len(g),
                    "von_baeyer": False,
                    "unit": [comps[i] for i in g],
                    "internal": sorted(s for s in internal if spiro[s][0] in g),
                }
            )
        for i in g:
            index_of[i] = len(merged) - 1
    outer = {}
    for s, members in spiro.items():
        if s in internal:
            continue
        mapped = sorted({index_of[m] for m in members})
        if len(mapped) != len(members):
            return None
        outer[s] = mapped
    return merged, outer


def _structure(mol, atoms=None):
    raw = _raw_structure(mol, atoms)
    if raw is None or not _connected(raw[1], len(raw[0])):
        return None
    comps, spiro = raw
    if _shape(spiro, len(comps))[0] is None:
        merged = _merge_units(comps, spiro)
        if merged is None:
            return None
        comps, spiro = merged
        if not _connected(spiro, len(comps)):
            return None
    neighbours = defaultdict(set)
    for members in spiro.values():
        for i in members:
            neighbours[i].update(m for m in members if m != i)
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


def has_spiro_union_shape(mol, atoms=None) -> bool:
    structure = _structure(mol, atoms)
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


def _unit_component(mol, comp):
    members = comp["unit"]
    if any(m["rings"] != 1 for m in members):
        raise UnsupportedStructure("a spiro unit of rings with several ring systems is not supported yet")
    for bond_idx in comp["bonds"]:
        if mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("unsaturated monocyclic spiro units are not supported yet")
    internal = comp["internal"]
    atom_rings = [m["atoms"] for m in members]
    if len(internal) != len(members) - 1:
        raise UnsupportedStructure("this spiro unit is not an unbranched chain of rings")
    count = len(members)
    holders = {s: [i for i, r in enumerate(atom_rings) if s in r] for s in internal}
    if any(len(h) != 2 for h in holders.values()):
        raise UnsupportedStructure("this spiro unit is not an unbranched chain of rings")
    degree = defaultdict(list)
    for s, (i, j) in holders.items():
        degree[i].append(j)
        degree[j].append(i)
    ends = [i for i in range(count) if len(degree[i]) == 1]
    if len(ends) != 2 or any(len(degree[i]) > 2 for i in range(count)):
        raise UnsupportedStructure("this spiro unit is not an unbranched chain of rings")
    order = [ends[0]]
    while len(order) < count:
        order.append(next(n for n in degree[order[-1]] if n not in order))
    chain_atoms = tuple(next(s for s, h in holders.items() if set(h) == {x, y}) for x, y in zip(order, order[1:]))
    graph = adjacency(mol)
    hetero = [a for a in comp["atoms"] if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    prefix = numerical_term(len(chain_atoms)) + "spiro" if len(chain_atoms) > 1 else "spiro"
    best, candidates = None, []
    for ring_order, spiros in _chain_direction_candidates(tuple(order), chain_atoms, atom_rings):
        options = _arc_choice_options(graph, atom_rings, ring_order, spiros)
        first_rest = atom_rings[ring_order[0]] - {spiros[0]}
        last_rest = atom_rings[ring_order[-1]] - {spiros[-1]}
        for start in (a for a in graph[spiros[0]] if a in first_rest):
            for end in (a for a in graph[spiros[-1]] if a in last_rest):
                for choices in product(*options):
                    seq, descriptor, superscripts = _build_sequence(graph, atom_rings, ring_order, spiros, start, end, choices)
                    if count == 2:
                        superscripts[-1] = None
                    locants = {a: str(i + 1) for i, a in enumerate(seq)}
                    descriptor_text = ".".join(
                        str(n) if sup is None else f"{n}^{locants[sup]}" for n, sup in zip(descriptor, superscripts)
                    )
                    parent = f"{prefix}[{descriptor_text}]{alkane_name(len(seq))}"
                    hetero_locants = sorted(int(locants[a]) for a in hetero)
                    key = (sorted(int(locants[x]) for x in spiros), tuple(descriptor), hetero_locants)
                    candidates.append((key, locants, parent))
    best = min(k for k, _, _ in candidates)
    kept = [(loc, parent) for k, loc, parent in candidates if k == best]
    parents = {parent for _, parent in kept}
    if len(parents) != 1:
        raise UnsupportedStructure("this spiro unit has no unique descriptor")
    (parent,) = parents

    def text(numbering):
        return (_replacement_prefixes(mol, hetero, numbering) if hetero else "") + parent

    return {
        "name": parent,
        "numberings": [(loc, ()) for loc, _ in kept],
        "replacement": False,
        "fused": False,
        "unit": True,
        "text": text,
    }


def _blocks(layout):
    group_of = {c: gi for gi, group in enumerate(layout["groups"]) for c in group}
    flat_index = {c: i for i, c in enumerate(layout["flat"])}
    blocks = defaultdict(list)
    for s, x, y in layout["links"]:
        early, late = sorted((x, y), key=lambda c: (group_of[c], flat_index[c]))
        blocks[max(group_of[early], group_of[late])].append((flat_index[late], flat_index[early], s, early, late))
    return {gi: [entry[2:] for entry in sorted(entries)] for gi, entries in blocks.items()}




def _analyse(mol, atoms=None):
    structure = _structure(mol, atoms)
    if structure is None:
        raise UnsupportedStructure("this is not a spiro union of rings")
    comps, spiro, _ = structure
    scope = set().union(*(c["atoms"] for c in comps))
    inside = [a for a in mol.GetAtoms() if atoms is None or a.GetIdx() in scope]
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in inside):
        raise UnsupportedStructure("isotopes, radicals and stereodescriptors of a spiro union are not supported yet")
    if any(
        b.GetStereo() != Chem.BondStereo.STEREONONE
        for b in mol.GetBonds()
        if atoms is None or (b.GetBeginAtomIdx() in scope and b.GetEndAtomIdx() in scope)
    ):
        raise UnsupportedStructure("double-bond stereo of a spiro union is not supported yet")
    info = {s: _spiro_atom(mol, s) for s in spiro}
    if any(a.GetFormalCharge() and a.GetIdx() not in spiro for a in inside):
        raise UnsupportedStructure("a charge away from the spiro atoms is not supported yet")
    cationic = [s for s, (_, charged) in info.items() if charged]
    if len(cationic) > 1:
        raise UnsupportedStructure("several cationic spiro atoms are not supported yet")
    if atoms is None and not cationic:
        for a in mol.GetAtoms():
            if a.GetIdx() not in scope and a.GetAtomicNum() != 6 and a.GetAtomicNum() not in _HALOGENS:
                raise UnsupportedStructure("a spiro union with a principal characteristic group is not supported yet")

    forced = set()
    for s, members in spiro.items():
        if mol.GetAtomWithIdx(s).GetAtomicNum() != 6 and not info[s][0]:
            if not any(comps[m]["von_baeyer"] for m in members):
                raise UnsupportedStructure("a standard-valence heteroatom at the spiro atom is not supported without a von Baeyer component")
            forced.update(members)
    named = [_unit_component(mol, c) if "unit" in c else _component(mol, c, i in forced) for i, c in enumerate(comps)]
    if any("unit" in c for c in comps) and any(
        mol.GetAtomWithIdx(s).GetAtomicNum() != 6 and info[s][0] and any("unit" in comps[m] for m in members)
        for s, members in spiro.items()
    ):
        raise UnsupportedStructure("a nonstandard spiro heteroatom of a spiro unit is not supported yet")

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
        if "text" in n:
            n["display"] = n["text"](n["numberings"][0][0])
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
    return SimpleNamespace(
        comps=comps, spiro=spiro, info=info, cationic=cationic, named=named, inner=inner, names=names, spiro_in=spiro_in,
        capable=capable, polycyclic_bonds=polycyclic_bonds, monocyclic_bonds=monocyclic_bonds,
        replacement_atoms=replacement_atoms, layouts=layouts, projections=projections, scope=scope,
    )


def _solutions(ctx, mol):
    indices = list(range(len(ctx.comps)))
    nonstandard = any(lam for lam, _ in ctx.info.values())
    for layout in ctx.layouts:
        flat = layout["flat"]
        prime_of = {c: k for k, c in enumerate(flat)}
        blocks = _blocks(layout)
        centre = layout.get("centre")
        terminal_atoms = [s for s, _, _ in layout["links"]] if centre is not None else []
        shortlist, shortlist_key = [], None
        for combo in product(*(list(ctx.projections[i]) for i in indices)):
            base = {}
            for i, proj in zip(indices, combo):
                for s, loc in zip(ctx.spiro_in[i], proj):
                    base[(s, i)] = loc
            locant = {k: _primed(loc, prime_of[k[1]]) for k, loc in base.items()}
            central = tuple(_lk(base[(s, centre)]) for s in terminal_atoms) if centre is not None and not nonstandard else ()
            citation = tuple(
                _lk(locant[(s, c)]) for gi in sorted(blocks) for s, early, late in blocks[gi] for c in (early, late)
            )
            key = (central, sorted(_lk(l) for l in locant.values()), citation)
            if shortlist_key is None or key < shortlist_key:
                shortlist, shortlist_key = [combo], key
            elif key == shortlist_key:
                shortlist.append(combo)
        for combo in shortlist:
            tables = [ctx.projections[i][proj] for i, proj in zip(indices, combo)]
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
                ene = [_primed(str(x), prime_of[i]) for i, (_, ene_i) in zip(indices, assignment) for x in ene_i]
                yield SimpleNamespace(
                    layout=layout, blocks=blocks, key=shortlist_key, assignment=assignment, locant_of=locant_of,
                    prime_of=prime_of, ene=ene,
                )


def _hetero_key(ctx, mol, locant_of):
    hetero = sorted(
        ctx.replacement_atoms, key=lambda a: (HETERO_RANK[mol.GetAtomWithIdx(a).GetSymbol()], _lk(locant_of[a]))
    )
    return (
        sorted(_lk(locant_of[a]) for a in hetero),
        [_lk(locant_of[a]) for a in hetero],
        [-(_bonding_number(mol.GetAtomWithIdx(a)) or 0) for a in sorted(hetero, key=lambda a: _lk(locant_of[a]))],
    )


def _ending_text(ene, ending):
    if len(ene) > 1:
        raise UnsupportedStructure("several double bonds of a spiro component are not supported yet")
    if ene and ending and ending[0] == "ylium":
        raise UnsupportedStructure("a double bond of a spiro cation is not supported yet")
    text, initial = "", ""
    if ene:
        text, initial = (f"-{ene[0]}-en" if ending else f"-{ene[0]}-ene"), "e"
    if ending:
        if ending[0] == "ylium":
            text += f"-{ending[1]}-ylium"
            initial = initial or "y"
        else:
            _, word, locants, added = ending
            added_text = f"({','.join(f'{a}H' for a in added)})" if added else ""
            text += f"-{','.join(map(str, locants))}{added_text}-{word}"
            initial = initial or word[0]
    return text, initial


def _compose(ctx, mol, sol, indicated, hydro, prefixes, ending):
    comps, spiro, info, named, inner = ctx.comps, ctx.spiro, ctx.info, ctx.named, ctx.inner
    layout, assignment, locant_of, blocks, prime_of = sol.layout, sol.assignment, sol.locant_of, sol.blocks, sol.prime_of
    flat, groups = layout["flat"], layout["groups"]

    def locant_in(s, c):
        return _primed(assignment[c][0][s], prime_of[c])

    def shown(c):
        text = named[c]["text"](assignment[c][0]) if "text" in named[c] else named[c]["name"]
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

    lam_front = ",".join(
        f"{lowest}λ{info[s][0]}"
        for lowest, s in sorted(
            ((min((locant_in(s, m) for m in spiro[s]), key=_lk), s) for s in spiro if info[s][0]),
            key=lambda pair: _lk(pair[0]),
        )
    )
    count = multiplying_prefix(len(spiro)) if len(spiro) > 1 else ""
    opening, closing = ("{", "}") if any("unit" in c for c in comps) else ("[", "]")
    kind = layout["kind"]
    identical_pair = kind == "chain" and len(comps) == 2 and len(set(ctx.names.values())) == 1
    identical_chain = kind == "chain" and len(comps) == 3 and len(set(ctx.names.values())) == 1
    if kind == "hub_ter":
        hub = layout["hub"]
        locants = sorted((locant_in(hub, c) for c in flat), key=_lk)
        lam = info[hub][0]
        last = shown(flat[0])
        body = f"{','.join([f'{locants[0]}λ{lam}' if lam else locants[0], *locants[1:]])}-spiroter{opening}{last}{closing}"
    elif identical_pair:
        first, second = flat
        spiro_atom = next(iter(spiro))
        lam = info[spiro_atom][0]
        last = shown(first)
        body = f"{locant_in(spiro_atom, first)}{f'λ{lam}' if lam else ''},{locant_in(spiro_atom, second)}-spirobi{opening}{last}{closing}"
    elif identical_chain:
        if lam_front:
            raise UnsupportedStructure("a nonstandard spiro atom in three identical components is not supported yet")
        pairs = ":".join(pair_text(gi) for gi in sorted(blocks))
        last = shown(flat[0])
        body = f"{pairs}-{count}spiroter{opening}{last}{closing}"
    else:
        parts = []
        for gi, group in enumerate(groups):
            if gi:
                parts.append(pair_text(gi))
            text = group_text(group)
            if gi and "unit" in comps[group[0]] and text[0].isdigit():
                text = re.sub(r"^([^-]+)-", r"[\1]", text)
            parts.append(text)
        last = parts[-1]
        if kind == "hub_different":
            parts[2] = f"({parts[2]})"
        body = f"{count}spiro{opening}{'-'.join(parts)}{closing}"
        if lam_front:
            body = f"{lam_front}-{body}"
    replacement_prefix = _replacement_prefixes(mol, ctx.replacement_atoms, locant_of) if ctx.replacement_atoms else ""
    if replacement_prefix:
        body = replacement_prefix + ("-" if body[0].isdigit() else "") + body
    ending_text, initial = _ending_text(sol.ene, ending)
    if ending_text:
        if initial in "aeiouy" and last.endswith("e") and not last.endswith(")"):
            body = body[: -len(last) - 1] + last[:-1] + closing
        body += ending_text

    out = prefixes
    for token in (hydro, indicated, body):
        if not token:
            continue
        joined = out and (token[0].isdigit() or token[0] in "[(" or out.endswith("H"))
        out += ("-" if joined else "") + token
    return out


def name_spiro_union(mol) -> str:
    from ._ylium_ring import _prefixes

    ctx = _analyse(mol)
    graph = adjacency(mol)
    ring_atoms = ctx.scope
    saturated, choices = _hydrogen_choices(mol, ctx.capable, ctx.polycyclic_bonds, ctx.monocyclic_bonds)
    best = None
    for sol in _solutions(ctx, mol):
        locant_of = sol.locant_of
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
        unsaturation = sorted(choice[1] + [_lk(x) for x in sol.ene])
        key = (sol.key, _hetero_key(ctx, mol, locant_of), choice[0], unsaturation, locant_set, citation_prefix)
        if best is None or key < best[0]:
            best = (key, sol, choice[2], grouped)
    if best is None:
        raise UnsupportedStructure("this spiro union has no supported numbering")
    _, sol, picked, grouped = best
    locant_of = sol.locant_of
    indicated = ",".join(f"{locant_of[a]}H" for a in sorted(picked, key=lambda a: _lk(locant_of[a])))
    hydro_atoms = sorted(saturated - picked, key=lambda a: _lk(locant_of[a]))
    hydro = f"{','.join(locant_of[a] for a in hydro_atoms)}-{multiplied_word(len(hydro_atoms), 'hydro')}" if hydro_atoms else ""
    prefixes = format_substituent_prefixes(grouped) if grouped else ""
    ending = None
    if ctx.cationic:
        s = ctx.cationic[0]
        ending = ("ylium", min((_primed(sol.assignment[m][0][s], sol.prime_of[m]) for m in ctx.spiro[s]), key=_lk))
    return _compose(ctx, mol, sol, indicated, hydro, prefixes, ending)


def spiro_union_numberings(mol, graph, skeleton_atoms):
    """Numberings of a spiro union of rings, with its name text for a suffix or free valence, for the engine that
    names ring systems bearing principal characteristic groups."""
    from ._ring_diyl_numbering import SUFFIX_ATOMS, Numbering, _exocyclic_oxo, _locs, _split_hydrogen, _yl

    atoms = set(skeleton_atoms)
    ctx = _analyse(mol, atoms)
    if ctx.cationic:
        raise UnsupportedStructure("a cationic spiro union is not named through its suffix groups")
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    doubled = set()
    for bond in kekule.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a not in ctx.scope or b not in ctx.scope:
            continue
        if bond.GetBondTypeAsDouble() == 2.0:
            if bond.GetIdx() in ctx.monocyclic_bonds:
                continue
            if a not in ctx.capable or b not in ctx.capable or bond.GetIdx() not in ctx.polycyclic_bonds:
                raise UnsupportedStructure("a double bond outside the polycyclic spiro components is not supported yet")
            doubled |= {a, b}
        elif bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("this bond order is not supported in a spiro system")
    saturated = ctx.capable - doubled
    skeleton = nx.Graph()
    skeleton.add_nodes_from(ctx.capable)
    adj = {a: set() for a in ctx.capable}
    for bond_idx in ctx.polycyclic_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in ctx.capable and b in ctx.capable:
            skeleton.add_edge(a, b)
            adj[a].add(b)
            adj[b].add(a)
    ih_count = len(skeleton) - 2 * len(nx.max_weight_matching(skeleton, maxcardinality=True))
    oxo_all = _exocyclic_oxo(mol, atoms) & ctx.capable
    suffix_atoms = SUFFIX_ATOMS.get() & atoms
    oxo_suffix = oxo_all & suffix_atoms
    accommodated = oxo_suffix | {a for a in suffix_atoms & saturated if mol.GetAtomWithIdx(a).GetAtomicNum() != 6}

    out = []
    for sol in _solutions(ctx, mol):
        position_of = {a: SpiroLocant(text) for a, text in sol.locant_of.items()}
        split = _split_hydrogen(position_of, adj, set(ctx.capable), set(saturated), oxo_all, accommodated, ih_count)
        if split is None:
            continue
        ih, added, hydro = split

        def text(locants, valence, substituted=frozenset(), suffix="yl", sol=sol, ih=ih, added=added, hydro=hydro):
            if suffix == "":
                ending = None
            elif suffix == "carboxylate":
                raise UnsupportedStructure("a carboxylate suffix on a spiro union is not supported yet")
            else:
                word = _yl(valence) if suffix == "yl" else multiplied_word(valence, suffix)
                ending = ("suffix", word, sorted(locants), added)
            hydro_text = f"{_locs(hydro)}-{multiplied_word(len(hydro), 'hydro')}" if hydro else ""
            indicated = ",".join(f"{p}H" for p in ih)
            return _compose(ctx, mol, sol, indicated, hydro_text, "", ending)

        numbering = Numbering(
            position_of,
            text,
            pre_key=(sol.key, _hetero_key(ctx, mol, sol.locant_of), ih),
            unsat_key=(tuple(sorted(SpiroLocant(x) for x in sol.ene)), added, hydro),
            ih=ih,
        )
        numbering.hydro, numbering.added = tuple(hydro), tuple(added)
        numbering.hydro_positions, numbering.added_positions = tuple(hydro), tuple(added)
        out.append(numbering)
    if not out:
        raise UnsupportedStructure("no numbering of this spiro union fits its hydro/indicated hydrogen")
    return out
