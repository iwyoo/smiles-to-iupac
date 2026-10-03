"""Metallocenes as parent hydrides (P-69.2.7): the retained -ocene names of
Fe/Ru/Os/Ni/Cr/Co/V sandwiches, substituted by prefixes and suffixes with the
ring atoms numbered 1-5 and 1'-5' ('1,1'-dimethylferrocene', 'ferrocen-1-ol',
'ferrocene-1,1'-dicarboxylic acid'), or cited as a substituent group
('1'-methylferrocen-1-yl') on a parent named by the other modules.
"""

import itertools

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._numerals import alkane_name, multiplying_prefix, numerical_term
from ._phosphanyl_group import PREFIX_PROP
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
    elif metal.GetDegree() == 0 and len(rings) == 2:
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


def _component_kind(mol, graph, atoms):
    """Most senior characteristic-group class (suffix candidate) in a
    substituent component, or None."""
    best = None
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
        if kind and (best is None or _ORDER.index(kind) < _ORDER.index(best)):
            best = kind
    return best


def _linker_name(graph, atoms, ends):
    if len(atoms) < 1 or len(ends) != 2:
        raise UnsupportedStructure("this linker between metallocenes is not supported here")
    inner = {a: [v for v in graph[a] if v in atoms] for a in atoms}
    if any(len(v) > 2 for v in inner.values()) or (len(atoms) > 1 and any(len(inner[e]) != 1 for e in ends)):
        raise UnsupportedStructure("this linker between metallocenes is not supported here")
    if sum(len(v) for v in inner.values()) != 2 * (len(atoms) - 1):
        raise UnsupportedStructure("a cyclic linker between metallocenes is not supported here")
    return f"{alkane_name(len(atoms))}-1,{len(atoms)}-diyl" if len(atoms) > 1 else "methylene"


def _name_phane(mol, units, graph, parent):
    """Ocenes joined 1,1'- in one ring through CH2 nodes (P-52.2.5):
    '1,3(1,1')-diferrocenacyclotetraphane'."""
    unit_atoms = [{m} | set(r[0]) | set(r[1]) for m, r in units]
    every = set().union(*unit_atoms)
    owner = {a: k for k, atoms in enumerate(unit_atoms) for a in atoms}
    roots = [(a, x) for m, rs in units for r in rs for a in r for x in graph[a] if x not in every and x != m]
    if len(roots) != 2 * len(units):
        return None
    ring_of = {a: ri for _, rs in units for ri, r in enumerate(rs) for a in r}
    comps = _components(graph, every, roots)
    links = {}
    for a, x in roots:
        comp = comps[x]
        if any(
            mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).GetIsAromatic()
            or mol.GetAtomWithIdx(i).GetTotalNumHs() != 2 or mol.GetAtomWithIdx(i).GetDegree() != 2
            for i in comp
        ):
            return None
        links.setdefault(frozenset(comp), []).append((owner[a], ring_of[a]))
    if any(len(e) != 2 or e[0][0] == e[1][0] for e in links.values()):
        return None
    exits = {}
    for comp, ends in links.items():
        for (u, r), (v, rv) in (ends, ends[::-1]):
            exits[(u, r)] = (v, rv, len(comp))
    cycle, unit, ring = [], 0, 0
    while True:
        v, rv, length = exits[(unit, ring)]
        cycle.append((unit, length))
        unit, ring = v, 1 - rv
        if unit == 0:
            break
    if len(cycle) != len(units) or ring != 0:
        return None
    total = sum(1 + n for _, n in cycle)
    lengths = [n for _, n in cycle]
    best = None
    for seq in (lengths, [lengths[(-k - 2) % len(lengths)] for k in range(len(lengths))]):
        for start in range(len(seq)):
            rotated = seq[start:] + seq[:start]
            positions, pos = [], 1
            for n in rotated:
                positions.append(pos)
                pos += 1 + n
            if best is None or positions < best:
                best = positions
    count = len(units)
    mult = multiplying_prefix(count)
    stem = parent[:-1] + "a"
    return f"{','.join(map(str, best))}(1,1{_PRIME})-{mult}{stem}cyclo{numerical_term(total)}phane"


def _name_paired_ocenes(mol, units, graph):
    parents = {OCENES[mol.GetAtomWithIdx(m).GetSymbol()] for m, _ in units}
    if len(parents) != 1 or any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("different metallocenes joined together are not supported here")
    parent = parents.pop()
    phane = _name_phane(mol, units, graph, parent)
    if phane:
        return phane
    if len(units) != 2:
        raise UnsupportedStructure("this arrangement of several metallocenes is not supported here")
    unit_atoms = [{m} | set(r[0]) | set(r[1]) for m, r in units]
    every = unit_atoms[0] | unit_atoms[1]
    roots = [(a, x, k) for k, (m, rs) in enumerate(units) for r in rs for a in r for x in graph[a] if x not in every and x != m]
    comps = _components(graph, every, [(a, x) for a, x, _ in roots])
    link_roots = [(a, x, k) for a, x, k in roots if any(v in unit_atoms[1 - k] for c in comps[x] for v in graph[c])]
    if len(link_roots) != 2 or {k for _, _, k in link_roots} != {0, 1}:
        raise UnsupportedStructure("this link between metallocenes is not supported here")
    linker = comps[link_roots[0][1]]
    if comps[link_roots[1][1]] != linker or any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).GetIsAromatic() for i in linker):
        raise UnsupportedStructure("this link between metallocenes is not supported here")
    link_name = _linker_name(graph, linker, [x for _, x, _ in link_roots])
    namer = PrefixNamer(mol, graph)
    results = []
    for k, (m, rs) in enumerate(units):
        unit_roots = [(a, x) for a, x, kk in roots if kk == k]
        named = {}
        for a, x in unit_roots:
            if x in linker:
                named[(a, x)] = (link_name, True)
            elif _component_kind(mol, graph, comps[x]) or _suffix_kind(mol, graph, a, x):
                raise UnsupportedStructure("characteristic groups on linked metallocenes are not supported here")
            else:
                named[(a, x)] = namer.name(x, a)
        best = None
        for positions in _numberings(graph, rs):
            grouped = _grouped(positions, named)
            key = (sorted(positions[a] for a, _ in named), [grouped[n]["locants"] for n in sorted(grouped)])
            if best is None or key < best[0]:
                best = (key, grouped)
        grouped = best[1]
        link_locant = grouped.pop(link_name)["locants"][0]
        results.append((link_locant, format_substituent_prefixes(grouped)))
    if results[0][1] != results[1][1]:
        raise UnsupportedStructure("different metallocene units joined together are not supported here")
    first, second = results[0][0], results[1][0]
    second += _PRIME * 2
    prefixes = results[0][1]
    multiplier = f"bis({prefixes}{parent})" if prefixes else f"di{parent}"
    return f"{first},{second}-({link_name}){multiplier}"


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
    if len(roots) == 2 and roots[0][1] != roots[1][1] and comps[roots[0][1]] & comps[roots[1][1]]:
        return _as_bridge(mol, unit, roots, comps, parent)
    for ring_atom, root in roots:
        for other_atom, other_root in roots:
            if other_root != root and comps[root] & comps[other_root]:
                raise UnsupportedStructure("a bridge between the rings of a metallocene is not supported here")
        if any(v in unit and v != ring_atom for a in comps[root] for v in graph[a]):
            raise UnsupportedStructure("a bridge between the rings of a metallocene is not supported here")
    suffix = {(a, x): _suffix_kind(mol, graph, a, x) for a, x in roots}
    comp_kind = {(a, x): _component_kind(mol, graph, comps[x]) for a, x in roots}
    ring_best = min((k for k in suffix.values() if k), key=_SUFFIX_ORDER.index, default=None)
    ring_rank = _ORDER.index(ring_best) if ring_best else len(_ORDER)
    chain_roots = [(a, x) for a, x in roots if not suffix[(a, x)] and comp_kind[(a, x)]]
    chain_best = min((comp_kind[r] for r in chain_roots), key=_ORDER.index, default=None)
    chain_rank = _ORDER.index(chain_best) if chain_best else len(_ORDER)
    namer = PrefixNamer(mol, graph)

    if chain_best and chain_rank < ring_rank:
        holders = [r for r in chain_roots if comp_kind[r] == chain_best]
        others = [r for r in roots if r not in holders]
        named = {}
        for a, x in others:
            kind = suffix[(a, x)]
            named[(a, x)] = (_PREFIX_WORD[kind], False) if kind else namer.name(x, a)
        attempt = None
        sp3 = all(
            mol.GetAtomWithIdx(x).GetAtomicNum() == 6
            and not mol.GetAtomWithIdx(x).GetIsAromatic()
            and all(b.GetBondTypeAsDouble() == 1.0 or _is_acyl_oxo(mol, x, b) for b in mol.GetAtomWithIdx(x).GetBonds())
            for _, x in holders
        )
        if not sp3:
            attempt = None
        elif len(holders) == 1:
            attempt = lambda: _as_substituent_group(mol, graph, rings, unit, named, holders[0], parent)
        elif len(holders) == 2 and not named:
            attempt = lambda: _as_multiplied_parent(mol, graph, rings, unit, comps, holders, parent)
        if attempt is not None:
            try:
                return attempt()
            except UnsupportedStructure:
                pass
        for a, x in holders:
            named[(a, x)] = namer.name(x, a)
        return _as_parent(mol, graph, rings, roots, suffix, None, named, parent)
    if chain_best and chain_rank == ring_rank:
        raise UnsupportedStructure("equally senior groups on the metallocene and on a substituent are not supported here")
    principal = ring_best
    named = {}
    for a, x in roots:
        kind = suffix[(a, x)]
        if kind and kind == principal:
            continue
        named[(a, x)] = (_PREFIX_WORD[kind], False) if kind else namer.name(x, a)
    return _as_parent(mol, graph, rings, roots, suffix, principal, named, parent)


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


def _contract(mol, graph, unit, named, main_root, group_name):
    remove = set(unit)
    for (a, x) in named:
        remove |= _group_atoms(graph, x, unit)
    rw = Chem.RWMol(mol)
    marker = rw.AddAtom(Chem.Atom(53))
    rw.AddBond(main_root, marker, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(marker).SetProp(PREFIX_PROP, enclose(group_name))
    for idx in sorted(remove, reverse=True):
        rw.RemoveAtom(idx)
    contracted = rw.GetMol()
    Chem.SanitizeMol(contracted)
    from .core import _name_mol

    name = _name_mol(contracted)
    if "iodo" in name or "iodide" in name or group_name not in name:
        raise UnsupportedStructure("the metallocene group could not be cited as a prefix here")
    return name


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


def _as_multiplied_parent(mol, graph, rings, unit, comps, holders, parent):
    """Two identical characteristic-group components joined by the metallocene
    ('1,1'-(ferrocene-1,1'-diyl)di(ethan-1-one)', P-15.3)."""
    (a1, x1), (a2, x2) = holders
    keys = [_component_smiles(mol, comps[x] | {a}, a) for a, x in holders]
    if len(comps[x1]) != len(comps[x2]) or keys[0] != keys[1]:
        raise UnsupportedStructure("different substituents on the two rings of a metallocene are not supported here")
    acyl = any(_is_acyl_oxo(mol, x1, b) for b in mol.GetAtomWithIdx(x1).GetBonds())
    rw = Chem.RWMol(mol)
    remove = set(unit) | comps[x2]
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
    inner = f"({stripped})" if acyl or any(ch.isdigit() for ch in stripped) else stripped
    ring_locs = sorted([(1, 0), (1, 1)])
    divalent = f"{parent}-{_fmt(ring_locs[0])},{_fmt(ring_locs[1])}-diyl"
    locs = f"{locant},{locant}{_PRIME}" if locant else ""
    head = f"{locs}-" if locs else ""
    return f"{head}({divalent})di{inner}"


def _as_bridge(mol, unit, roots, comps, parent):
    """A chain joining both rings, cited as a divalent metallocene group
    ('3,5-(ferrocene-1,1'-diyl)pentanoic acid')."""
    if roots[0][0] in roots[1][0:1] or len({a for a, _ in roots}) != 2:
        raise UnsupportedStructure("this metallocene bridge is not supported here")
    token = "QQQ"
    rw = Chem.RWMol(mol)
    for _, x in roots:
        marker = rw.AddAtom(Chem.Atom(53))
        rw.AddBond(x, marker, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(marker).SetProp(PREFIX_PROP, token)
    for idx in sorted(unit, reverse=True):
        rw.RemoveAtom(idx)
    contracted = rw.GetMol()
    Chem.SanitizeMol(contracted)
    from .core import _name_mol
    import re

    name = _name_mol(contracted)
    match = re.match(r"(\d+),(\d+)-(?:di|bis\()" + token + r"\)?", name)
    if not match or name.count(token) != 1 or "iod" in name.replace(token, ""):
        raise UnsupportedStructure("the metallocene bridge could not be cited as a divalent group here")
    first, second = match.group(1), match.group(2)
    rest = name[match.end():].lstrip("-")
    return f"{first},{second}-({parent}-1,1{_PRIME}-diyl){rest}"


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
