"""Ring assemblies of identical monocycles named on their longest unbranched chain of three to six rings, the other
rings becoming substituents (P-28.3, P-28.6): [11,21:24,31-terphenyl]-14-ol, 25-phenyl-11,21:23,31-terphenyl,
[11,21:24,31-tercyclohexan]-14-yl."""

import re
from itertools import product

from ._common import multiplied_word, ring_cycle
from ._multiplicative import _bare_key
from ._multiplicative_ring import _SUFFIX_WORDS, _citation_key, monocycle_spec, numberings
from ._multiplicative_text import CompositeLocant
from ._substituents import format_substituent_prefixes, name_branch

_LATIN = {3: "ter", 4: "quater", 5: "quinque", 6: "sexi"}
_MAX_RINGS = 10
_MAX_FUSED_UNITS = 4


def _cite(locant):
    return f"{locant[0]}{locant[1]}"


def _longest_paths(mol, rings):
    """Longest simple paths through the tree of rings as (ring order, junction atom pairs); None when the rings
    are not a tree of single bonds or the longest path exceeds six rings."""
    owner = {a: i for i, ring in enumerate(rings) for a in ring}
    links = {i: [] for i in range(len(rings))}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in owner and b in owner and owner[a] != owner[b]:
            if bond.GetBondTypeAsDouble() != 1.0 or bond.GetIsAromatic():
                return None
            links[owner[a]].append((owner[b], a, b))
            links[owner[b]].append((owner[a], b, a))
    if sum(len(v) for v in links.values()) != 2 * (len(rings) - 1):
        return None
    if any(len({x[0] for x in v}) != len(v) for v in links.values()):
        return None
    found = []

    def walk(order, junctions):
        for nxt, a, b in links[order[-1]]:
            if nxt not in order:
                walk(order + [nxt], junctions + [(a, b)])
        found.append((order, junctions))

    for start in links:
        walk([start], [])
    longest = max(len(o) for o, _ in found)
    if longest > 6 or longest < 3:
        return None
    result, seen = [], set()
    for order, junctions in found:
        if len(order) == longest and tuple(order[::-1]) not in seen:
            seen.add(tuple(order))
            result.append((order, junctions))
    return result


def _unit_systems(mol, graph, frees, within):
    from ._polyfunctional import _arm_atoms
    from ._system_assembly import _systems

    systems = _systems(mol)
    if frees:
        arm = within if within is not None else _arm_atoms(graph, frees[0][0], frees[0][1])
        systems = [(rings, atoms) for rings, atoms in systems if set(atoms) <= arm]
    return systems


def chain_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, free=None, within=None):
    """(count, ((-count,), name, parts)) for the parent form, (name, True) for the substituent form `free`
    = (root, coming_from) or a list of such pairs (within the atom set `within`); None when `mol` is not such
    an assembly."""
    from ._polyfunctional import _require_mancude_system, _stereo_rank
    from ._ring_diyl_numbering import SUFFIX_ATOMS, system_numberings
    from ._system_assembly import _fully_hydro, _hydro, _lit_hydrogens, _stem

    frees = [] if free is None else [free] if isinstance(free[0], int) else list(free)
    systems = _unit_systems(mol, graph, frees, within)
    rings = [sorted(atoms) for _, atoms in systems]
    if not 3 <= len(rings) <= _MAX_RINGS:
        return None
    fused = [len(r) > 1 for r, _ in systems]
    if any(fused) and not all(fused):
        return None
    if not frees and principal is not None and not occurrences:
        return None
    spec = None
    if not any(fused):
        specs = [monocycle_spec(mol, r) for r in rings]
        if any(s is None or s.kind == "pyrrole" for s in specs):
            return None
        if len({(s.kind, len(r)) for s, r in zip(specs, rings)}) != 1:
            return None
        spec = specs[0]
    elif len(rings) > _MAX_FUSED_UNITS:
        return None
    if len({_bare_key(mol, set(r)) for r in rings}) != 1:
        return None
    paths = _longest_paths(mol, rings)
    if not paths:
        return None
    if any(fused):
        for _, atoms in systems:
            _require_mancude_system(mol, atoms)
    owned = set().union(*(o[2] for o in occurrences)) if occurrences else set()
    marked = [o[1] for o in occurrences] if not frees else [f[0] for f in frees]
    best = None
    for order, junctions in paths:
        chain_rings = [rings[i] for i in order]
        chain_systems = [systems[i] for i in order]
        atoms_all = {a for r in chain_rings for a in r}
        if any(a not in atoms_all for a in marked):
            continue
        roots = [
            (r, n.GetIdx())
            for r in atoms_all
            for n in mol.GetAtomWithIdx(r).GetNeighbors()
            if n.GetIdx() not in atoms_all
            and n.GetIdx() not in owned
            and (r, n.GetIdx()) not in frees
        ]
        entries = [(r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots]
        ring_graph = {
            a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in atoms_all] for a in atoms_all
        }

        def orientations(unit, ring_atoms, attach):
            if fused[0]:
                token = SUFFIX_ATOMS.set(frozenset(attach) | (set(marked) & set(ring_atoms)))
                try:
                    return [(n.position_of, n) for n in system_numberings(mol, graph, unit[0], unit[1])]
                finally:
                    SUFFIX_ATOMS.reset(token)
            if spec.hetero is not None:
                return [(n, None) for n in numberings(monocycle_spec(mol, ring_atoms))]
            cycle = ring_cycle(ring_graph, ring_atoms)
            found = []
            for start in attach:
                rotated = cycle[cycle.index(start):] + cycle[: cycle.index(start)]
                for sequence in (rotated, [rotated[0]] + rotated[:0:-1]):
                    found.append(({atom: i + 1 for i, atom in enumerate(sequence)}, None))
            return found

        for direction in (0, 1):
            chain = chain_rings if direction == 0 else chain_rings[::-1]
            units = chain_systems if direction == 0 else chain_systems[::-1]
            pairs = junctions if direction == 0 else [(b, a) for a, b in junctions[::-1]]
            attach = [set() for _ in chain]
            for i, (a, b) in enumerate(pairs):
                attach[i].add(a)
                attach[i + 1].add(b)
            candidates = [orientations(units[i], chain[i], attach[i]) for i in range(len(chain))]
            for combo in product(*candidates):
                locants = {atom: (i + 1, pos) for i, (numbering, _) in enumerate(combo) for atom, pos in numbering.items()}
                if any(a not in locants for a in marked) or any(a not in locants for pair in pairs for a in pair):
                    continue
                head = ((), (), ())
                if fused[0]:
                    infos = [info for _, info in combo]
                    stems = {_stem(info) for info in infos}
                    hydros = {_hydro(info) for info in infos}
                    if len(stems) != 1 or None in stems or len(hydros) != 1:
                        continue
                    lit = sorted(
                        (i + 1, p) for i, info in enumerate(infos) for p in _lit_hydrogens(mol, info)
                    )
                    head = (tuple(lit), (stems.pop(), hydros.pop(), _fully_hydro(infos[0])), ())
                junction_text = [(locants[a], locants[b]) for a, b in pairs]
                flat = [x for pair in junction_text for x in pair]
                key = (
                    tuple(sorted(flat)),
                    tuple(flat),
                    head[0],
                    tuple(sorted(locants[a] for a in marked)),
                    -len(entries),
                    tuple(sorted(locants[r] for r, _, _ in entries)),
                    _citation_key([(locants[r], name) for r, name, _ in entries]),
                    _stereo_rank(stereo, {a: CompositeLocant(*loc) for a, loc in locants.items()}, True),
                )
                if best is None or key < best[0]:
                    best = (key, locants, junction_text, entries, len(chain), head)
    if best is None:
        return None
    _, locants, junction_text, entries, count, head = best
    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(_cite(locants[r]))
    for info in grouped.values():
        info["locants"].sort(key=lambda text: (int(text[0]), float(text[1:].rstrip("abcdefghijklmnopqrstuvwxyz") or 0), text))
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    junction_str = ":".join(f"{_cite(a)},{_cite(b)}" for a, b in junction_text)
    ih_text = ""
    if fused[0]:
        stem, hydro, full = head[1]
        ring_word = stem
        if hydro:
            cited = [f"{ring}{p}" for ring in range(1, count + 1) for p in hydro]
            word = multiplied_word(len(cited), "hydro")
            ih_text = (word if full else f"{','.join(cited)}-{word}") + "-"
        if head[0]:
            ih_text += ",".join(f"{ring}{p}H" for ring, p in head[0]) + "-"
    else:
        ring_word = "phenyl" if spec.kind == "benzene" else spec.parent

    def parent(elide):
        word = ring_word[:-1] if elide and ring_word.endswith("e") else ring_word
        word = f"({word})" if re.search(r"[\d\[-]", word) else word
        return f"{junction_str}-{_LATIN[count]}{word}"

    def base(elide):
        return ih_text + parent(elide)

    def bracketed(elide):
        return f"{ih_text}[{parent(elide)}]"

    if frees:
        spots = ",".join(_cite(c) for c in sorted(locants[f[0]] for f in frees))
        core = f"{bracketed(True)}-{spots}-{multiplied_word(len(frees), 'yl')}"
        return (f"{prefix}-{core}" if prefix else core), True
    total = len(occurrences)
    if principal is None:
        core = base(False)
    else:
        from ._polyfunctional import _RING_SUFFIX

        word = multiplied_word(total, _SUFFIX_WORDS[_RING_SUFFIX[principal]])
        spots = ",".join(_cite(locants[o[1]]) for o in sorted(occurrences, key=lambda o: locants[o[1]]))
        core = f"{bracketed(word[0] in 'aeiouy')}-{spots}-{word}"
    name = f"{prefix}-{core}" if prefix else core
    return total, ((-total,), name, (None, None, None, 0, {a: CompositeLocant(*loc) for a, loc in locants.items()}, True))


def assembly_diyl(mol, graph, halogens, aromatic_atoms, atoms, frees):
    """Name of the polyvalent group formed by identical monocycles joined directly ([1,1'-biphenyl]-4,4'-diyl),
    `frees` = [(attached ring atom, outside neighbour)]; None when the atoms are not such an assembly."""
    inside = set(atoms)
    rings = [list(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= inside]
    from ._system_assembly import system_assembly

    orders = {f[0]: int(mol.GetBondBetweenAtoms(*f).GetBondTypeAsDouble()) for f in frees}
    mixed = len(set(orders.values())) > 1
    if mixed and len(rings) != 2:
        return None
    fused = system_assembly(mol, graph, halogens, aromatic_atoms, None, [], None, free=frees, within=inside)
    if fused is not None:
        return fused[0]
    if len(rings) >= 3:
        found = chain_assembly(mol, graph, halogens, aromatic_atoms, None, [], None, free=frees, within=inside)
        return found[0] if found else None
    if len(rings) != 2:
        return None
    from ._polyfunctional import _assembly_base, _assembly_numbering, _locant_order

    specs = [monocycle_spec(mol, r) for r in rings]
    if any(s is None or s.kind == "pyrrole" for s in specs) or specs[0].kind != specs[1].kind or len(rings[0]) != len(rings[1]):
        return None
    if _bare_key(mol, set(rings[0])) != _bare_key(mol, set(rings[1])):
        return None
    joins = [(a, b) for a in rings[0] for b in rings[1] if mol.GetBondBetweenAtoms(a, b) is not None]
    if len(joins) != 1 or mol.GetBondBetweenAtoms(*joins[0]).GetBondTypeAsDouble() != 1.0:
        return None
    roots = [
        (r, n.GetIdx())
        for r in inside
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in inside and (r, n.GetIdx()) not in frees
    ]
    entries = [(r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots]
    marked = [f[0] for f in frees]
    cite_marked = sorted(marked, key=lambda a: orders[a]) if mixed else None
    locants = _assembly_numbering(graph, rings, joins[0], marked, entries, specs, cite_marked)

    def cite(locant):
        return f"{locant[1]}{chr(39) * locant[0]}"

    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(cite(locants[r]))
    for info in grouped.values():
        info["locants"].sort(key=lambda text: (int(text.rstrip(chr(39))), text.count(chr(39))))
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    base = _assembly_base(specs, locants, joins[0], elide=True)
    if mixed:
        words = {1: "yl", 2: "ylidene", 3: "ylidyne"}
        if any(order not in words for order in orders.values()):
            return None
        pieces = []
        for order in sorted(set(orders.values())):
            group = [a for a in marked if orders[a] == order]
            spots = ",".join(cite(c) for c in sorted((locants[m] for m in group), key=_locant_order))
            pieces.append(f"{spots}-{multiplied_word(len(group), words[order])}")
        core = f"[{base}]-" + "-".join(pieces)
        return f"{prefix}-{core}" if prefix else core
    spots = ",".join(cite(c) for c in sorted((locants[m] for m in marked), key=_locant_order))
    core = f"[{base}]-{spots}-{multiplied_word(len(frees), 'yl')}"
    return f"{prefix}-{core}" if prefix else core
