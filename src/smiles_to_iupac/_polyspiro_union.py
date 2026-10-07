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
    if any(len(members) != 2 for members in spiro.values()) or len(spiro) != len(comps) - 1:
        return None
    neighbours = defaultdict(set)
    for i, j in spiro.values():
        neighbours[i].add(j)
        neighbours[j].add(i)
    seen, queue = {0}, [0]
    while queue:
        for n in neighbours[queue.pop()]:
            if n not in seen:
                seen.add(n)
                queue.append(n)
    if len(seen) != len(comps):
        return None
    return comps, spiro, neighbours


def has_polyspiro_union_shape(mol) -> bool:
    structure = _structure(mol)
    if structure is None:
        return False
    _, _, neighbours = structure
    degrees = sorted(len(n) for n in neighbours.values())
    return degrees[-1] <= 2 or (degrees[-2] == 1 and degrees[-1] >= 3)


def _chain_layouts(neighbours, names):
    start = next(c for c in neighbours if len(neighbours[c]) == 1)
    path, previous = [start], None
    while len(path) < len(neighbours):
        following = next(n for n in neighbours[path[-1]] if n != previous)
        previous = path[-1]
        path.append(following)
    orientations = [path, path[::-1]]
    keyed = [([_name_key(names[c]) for c in o], o) for o in orientations]
    best = min(k for k, _ in keyed)
    chosen = [o for k, o in keyed if k == best]
    unique = []
    for o in chosen:
        if o not in unique:
            unique.append(o)
    return [{"kind": "chain", "flat": o, "groups": [[c] for c in o]} for o in unique]


def _classes(members, names):
    classes = defaultdict(list)
    for c in members:
        classes[names[c]].append(c)
    return sorted(classes.values(), key=lambda cls: _name_key(names[cls[0]]))


def _star_layouts(neighbours, names):
    centre = next(c for c in neighbours if len(neighbours[c]) >= 3)
    terminals = sorted(neighbours[centre])
    classes = _classes(terminals, names)
    layouts = []
    if len(classes) == 1:
        for order in permutations(terminals):
            layouts.append({"kind": "star", "flat": [centre, *order], "groups": [[centre], list(order)], "centre": centre})
        return layouts
    first, rest = classes[0], classes[1:]
    for first_order in permutations(first):
        for rest_orders in product(*(permutations(cls) for cls in rest)):
            flat = [first_order[0], centre, *first_order[1:], *(c for order in rest_orders for c in order)]
            groups = [list(first_order), [centre], *(list(order) for order in rest_orders)]
            layouts.append({"kind": "star", "flat": flat, "groups": groups, "centre": centre})
    return layouts


def _blocks(layout, spiro):
    group_of = {c: gi for gi, group in enumerate(layout["groups"]) for c in group}
    flat_index = {c: i for i, c in enumerate(layout["flat"])}
    blocks = defaultdict(list)
    for s, members in spiro.items():
        early, late = sorted(members, key=flat_index.get)
        blocks[max(group_of[early], group_of[late])].append((flat_index[late], s, early, late))
    return {gi: sorted(entries) for gi, entries in blocks.items()}


def name_polyspiro_union(mol) -> str:
    from ._ylium_ring import _prefixes

    comps, spiro, neighbours = _structure(mol)
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
            if lam - 4 > 1:
                raise UnsupportedStructure("a spiro atom with several multiple bonds is not supported yet")
            if lam - 4 == 1 and any(named[m]["fused"] for m in spiro[s]):
                capable.add(s)
    saturated, choices = _hydrogen_choices(mol, capable, polycyclic_bonds, monocyclic_bonds)

    degrees = {c: len(neighbours[c]) for c in neighbours}
    layouts = _chain_layouts(neighbours, names) if max(degrees.values()) <= 2 else _star_layouts(neighbours, names)

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
        blocks = _blocks(layout, spiro)
        centre = layout.get("centre")
        terminal_atoms = []
        if centre is not None:
            terminal_atoms = [
                s for c in flat if c != centre for s, members in spiro.items() if c in members
            ]
        indices = list(range(len(comps)))
        shortlist, shortlist_key = [], None
        for combo in product(*(list(projections[i]) for i in indices)):
            base = {}
            for i, proj in zip(indices, combo):
                for s, loc in zip(spiro_in[i], proj):
                    base[(s, i)] = loc
            locant = {k: _primed(loc, prime_of[k[1]]) for k, loc in base.items()}
            central = tuple(_lk(base[(s, centre)]) for s in terminal_atoms) if centre is not None else ()
            citation = tuple(
                _lk(locant[(s, c)]) for gi in sorted(blocks) for _, s, early, late in blocks[gi] for c in (early, late)
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
        return ":".join(f"{locant_in(s, early)},{locant_in(s, late)}" for _, s, early, late in blocks[gi])

    indicated = ",".join(f"{locant_of[a]}H" for a in sorted(picked, key=lambda a: _lk(locant_of[a])))
    hydro_atoms = sorted(saturated - picked, key=lambda a: _lk(locant_of[a]))
    hydro = f"{','.join(locant_of[a] for a in hydro_atoms)}-{multiplied_word(len(hydro_atoms), 'hydro')}" if hydro_atoms else ""
    prefixes = format_substituent_prefixes(grouped) if grouped else ""

    lam_front = ",".join(
        f"{lowest}λ{info[s][0]}"
        for lowest, s in sorted(
            (min((locant_in(s, m) for m in spiro[s]), key=_lk), s) for s in spiro if info[s][0]
        )
    )
    count = multiplying_prefix(len(spiro))
    identical_chain = layout["kind"] == "chain" and len(comps) == 3 and len({names[c] for c in range(3)}) == 1
    if identical_chain:
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
        body = f"{count}spiro[{'-'.join(parts)}]"
        if lam_front:
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
