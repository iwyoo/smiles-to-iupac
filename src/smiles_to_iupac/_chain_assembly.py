"""Ring assemblies of identical monocycles named on their longest unbranched chain of three to six rings, the other
rings becoming substituents (P-28.3, P-28.6): [11,21:24,31-terphenyl]-14-ol, 25-phenyl-11,21:23,31-terphenyl,
[11,21:24,31-tercyclohexan]-14-yl."""

import re
from itertools import product

from ._free_valence import SUFFIX_OF_ORDER, valence_word
from ._common import KEEP_ASSEMBLY_ENDING
from ._common import assembly_join, multiplied_word, ring_cycle, superscript_locant
from ._multiplicative import _bare_key
from ._multiplicative_ring import _SUFFIX_WORDS, _citation_key, numberings
from ._multiplicative_text import CompositeLocant
from ._ring_assembly_chain import _hetero_ring_alignments
from ._substituents import format_substituent_prefixes, name_branch

_LATIN = {3: "ter", 4: "quater", 5: "quinque", 6: "sexi"}
_MAX_RINGS = 10
_MAX_FUSED_UNITS = 4


def _cite(locant):
    return superscript_locant(*locant)


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


def _hydro_component_rings(mol, graph, rings):
    """The rings of an assembly that are the hydro forms of its mancude components, which then take hydro prefixes
    (P-54.3): benzene beside cyclohexane, cyclohexene or cyclohexadiene, pyridine beside piperidine; empty otherwise."""
    from ._ring_assembly_chain import ring_kinds_with_hydro

    kinds, hydro = ring_kinds_with_hydro(mol, graph, [tuple(r) for r in rings])
    if None in kinds or len(set(kinds)) != 1 or not hydro or len(hydro) == len(rings):
        return []
    return hydro


def _component_hydro(mol, locants, hydro_rings):
    """The sorted (ring number, position) of the atoms taking hydro prefixes: every position of a saturated or partly
    saturated six-membered ring, every one but the heteroatom of a five-membered ring (P-31.1.4.2.4)."""
    from ._ring_assembly_chain import hydrogenated_atoms

    found = []
    for ring in hydro_rings:
        if any(atom not in locants for atom in ring):
            continue
        for atom in hydrogenated_atoms(mol, ring):
            ring_number, position = locants[atom]
            if len(ring) == 5 and position == 1:
                continue
            found.append((ring_number, position))
    return tuple(sorted(found))


def _unit_systems(mol, graph, frees, within):
    from ._polyfunctional import _arm_atoms
    from ._system_assembly import _systems

    systems = _systems(mol)
    arm = within if within is not None else _arm_atoms(graph, frees[0][0], frees[0][1]) if frees else None
    if arm is not None:
        systems = [(rings, atoms) for rings, atoms in systems if set(atoms) <= arm]
    return systems


def chain_assembly(
    mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, free=None, within=None, place=None, n_names=()
):
    """(count, ((-count,), name, parts)) for the parent form, (name, True) for the substituent form `free`
    = (root, coming_from) or a list of such pairs (within the atom set `within`); None when `mol` is not such
    an assembly."""
    from ._diester_ring_diyl import _isotope_key
    from ._polyfunctional import _assembly_spec, _require_mancude_system, _stereo_rank
    from ._ring_diyl_numbering import LINK_ATOMS, SUFFIX_ATOMS, system_numberings
    from ._system_assembly import _fully_hydro, _hetero_mancude_ring, _hydro, _monocycle_numberings, _stem, unit_indicated_hydrogens

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
    hydro_rings = []
    if not any(fused):
        hydro_rings = _hydro_component_rings(mol, graph, rings)
        if hydro_rings:
            spec = _assembly_spec(mol, graph, next(r for r in rings if frozenset(r) not in hydro_rings))
            if spec is None or spec.kind == "pyrrole":
                return None
        else:
            specs = [_assembly_spec(mol, graph, r) for r in rings]
            if any(s is None or s.kind == "pyrrole" for s in specs):
                return None
            if len({(s.kind, len(r)) for s, r in zip(specs, rings)}) != 1:
                return None
            spec = specs[0]
    elif len(rings) > _MAX_FUSED_UNITS:
        return None
    from ._system_assembly import _skeleton_key

    if not hydro_rings and len({(_skeleton_key if any(fused) else _bare_key)(mol, set(r)) for r in rings}) != 1:
        return None
    paths = _longest_paths(mol, rings)
    if not paths:
        return None
    if any(fused):
        for _, atoms in systems:
            _require_mancude_system(mol, atoms)
    owned = set().union(*(o[2] for o in occurrences)) if occurrences else set()
    marked = [o[1] for o in occurrences] if not frees else [f[0] for f in frees]
    suffix_atoms = frozenset(o[1] for o in occurrences)
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
                links = LINK_ATOMS.set(frozenset(attach))
                try:
                    return [(n.position_of, n) for n in system_numberings(mol, graph, unit[0], unit[1])]
                finally:
                    LINK_ATOMS.reset(links)
                    SUFFIX_ATOMS.reset(token)
            if spec.hetero is not None:
                if frozenset(ring_atoms) in hydro_rings:
                    return [(n, None) for n in _hetero_ring_alignments(mol, graph, ring_atoms, spec.parent)]
                if suffix_atoms & set(ring_atoms) and _hetero_mancude_ring(mol, set(ring_atoms)):
                    token = SUFFIX_ATOMS.set(suffix_atoms & set(ring_atoms))
                    links = LINK_ATOMS.set(frozenset(attach))
                    try:
                        return [(n.position_of, n) for n in _monocycle_numberings(mol, graph, ring_atoms)]
                    finally:
                        LINK_ATOMS.reset(links)
                        SUFFIX_ATOMS.reset(token)
                return [(n, None) for n in numberings(_assembly_spec(mol, graph, ring_atoms))]
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
                infos = [info for _, info in combo]
                carries = [info is not None and bool(suffix_atoms & set(info.position_of)) for info in infos]
                if fused[0] or any(carries):
                    stem = None
                    if fused[0]:
                        stems = {_stem(info, carries[i]) for i, info in enumerate(infos)}
                        if len(stems) != 1 or None in stems:
                            continue
                        stem = stems.pop()
                    lit = sorted(
                        (i + 1, p)
                        for i, info in enumerate(infos)
                        if info is not None
                        for p in unit_indicated_hydrogens(mol, info, attach[i], carries[i], not fused[0])
                    )
                    units_hydro = [
                        (_hydro(info) if fused[0] else tuple(info.hydro_positions) if carries[i] else ()) if info is not None else ()
                        for i, info in enumerate(infos)
                    ]
                    hydro = tuple(sorted((i + 1, p) for i, found in enumerate(units_hydro) for p in found))
                    even = len(set(units_hydro)) == 1
                    full = even and all(info is not None and _fully_hydro(info) for info in infos)
                    added = tuple(sorted((i + 1, p) for i, info in enumerate(infos) if carries[i] for p in info.added_positions))
                    head = (tuple(lit), (stem, hydro, full), added)
                if hydro_rings:
                    own = head[1][1] if head[1] else ()
                    merged = tuple(sorted(set(own) | set(_component_hydro(mol, locants, hydro_rings))))
                    head = (head[0], (head[1][0] if head[1] else None, merged, False), head[2])
                junction_text = [(locants[a], locants[b]) for a, b in pairs]
                flat = [x for pair in junction_text for x in pair]
                key = (
                    tuple(sorted(flat)),
                    tuple(flat),
                    head[0],
                    tuple(sorted(locants[a] for a in marked)),
                    head[2],
                    head[1][1] if head[1] else (),
                    -len(entries),
                    tuple(sorted(locants[r] for r, _, _ in entries)),
                    _citation_key([(locants[r], name) for r, name, _ in entries]),
                    _stereo_rank(stereo, {a: CompositeLocant(*loc) for a, loc in locants.items()}, True),
                    _isotope_key(locants, locants),
                )
                if best is None or key < best[0]:
                    best = (key, locants, junction_text, entries, len(chain), head)
    if best is None:
        return None
    _, locants, junction_text, entries, count, head = best
    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locants[r])
    for info in grouped.values():
        info["locants"].sort(key=lambda loc: (loc[0], float(str(loc[1]).rstrip("abcdefghijklmnopqrstuvwxyz") or 0), str(loc[1])))
        info["locants"] = [_cite(loc) for loc in info["locants"]]
    if n_names:
        from ._polyfunctional import assembly_n_names

        grouped = assembly_n_names(grouped, n_names, {a: CompositeLocant(*loc) for a, loc in locants.items()}, len(occurrences))
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    junction_str = ":".join(f"{_cite(a)},{_cite(b)}" for a, b in junction_text)
    ih_text = ""
    if fused[0]:
        ring_word = head[1][0]
    else:
        ring_word = "phenyl" if spec.kind == "benzene" else spec.parent
    if head[1]:
        _, hydro, full = head[1]
        if hydro:
            cited = [superscript_locant(ring, p) for ring, p in hydro]
            word = multiplied_word(len(cited), "hydro")
            ih_text = (word if full else f"{','.join(cited)}-{word}") + "-"
    if head[0]:
        ih_text += ",".join(f"{superscript_locant(ring, p)}H" for ring, p in head[0]) + "-"

    def parent(elide):
        word = ring_word[:-1] if elide and not KEEP_ASSEMBLY_ENDING.get() and ring_word.endswith("e") else ring_word
        word = f"({word})" if re.search(r"[\d\[-]", word) else word
        return f"{junction_str}-{_LATIN[count]}{word}"

    def base(elide):
        return ih_text + parent(elide)

    def bracketed(elide):
        lead = ih_text[:-1] if ih_text.endswith("hydro-") else ih_text
        return f"{lead}[{parent(elide)}]"

    if frees:
        spots = ",".join(_cite(c) for c in sorted(locants[f[0]] for f in frees))
        core = f"{bracketed(True)}-{spots}-{multiplied_word(len(frees), 'yl')}"
        name = assembly_join(prefix, core)
        if place is not None:
            place.update(
                locants={a: CompositeLocant(*loc) for a, loc in locants.items()},
                start=len(name) - len(core),
                atoms=atoms_all,
            )
        return name, True
    total = len(occurrences)
    if principal is None:
        core = base(False)
    else:
        from ._polyfunctional import _RING_SUFFIX

        word = multiplied_word(total, _SUFFIX_WORDS[_RING_SUFFIX[principal]])
        spots = ",".join(_cite(locants[o[1]]) for o in sorted(occurrences, key=lambda o: locants[o[1]]))
        if head[2]:
            spots += "(" + ",".join(f"{superscript_locant(ring, p)}H" for ring, p in head[2]) + ")"
        core = f"{bracketed(word[0] in 'aeiouy')}-{spots}-{word}"
    name = assembly_join(prefix, core)
    return total, ((-total,), name, (None, None, None, 0, {a: CompositeLocant(*loc) for a, loc in locants.items()}, True, len(name) - len(core)))


def assembly_diyl(mol, graph, halogens, aromatic_atoms, atoms, frees, stereo=None):
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
        return None if stereo else fused[0]
    if len(rings) >= 3:
        found = chain_assembly(mol, graph, halogens, aromatic_atoms, None, [], None, free=frees, within=inside)
        if found is not None:
            return None if stereo else found[0]
    if not mixed:
        from ._vb_ring_assembly import vb_assembly

        von_baeyer = vb_assembly(mol, graph, halogens, aromatic_atoms, None, [], None, free=frees, within=inside)
        if von_baeyer is not None:
            return None if stereo else von_baeyer[0]
    if len(rings) != 2:
        return None
    from ._polyfunctional import (
        _assembly_base,
        _assembly_numbering,
        _assembly_stereo_prefix,
        _assembly_replacement,
        _assembly_spec,
        _locant_order,
        _replacement_junction_ok,
    )

    specs = [_assembly_spec(mol, graph, r, True) for r in rings]
    if any(s is None or s.kind == "pyrrole" for s in specs) or specs[0].kind != specs[1].kind or len(rings[0]) != len(rings[1]):
        return None
    if specs[0].name or specs[1].name:
        if specs[0].name != specs[1].name:
            return None
    elif not any(s.replaced for s in specs) and _bare_key(mol, set(rings[0])) != _bare_key(mol, set(rings[1])):
        return None
    joins = [(a, b) for a in rings[0] for b in rings[1] if mol.GetBondBetweenAtoms(a, b) is not None]
    if len(joins) != 1 or mol.GetBondBetweenAtoms(*joins[0]).GetBondTypeAsDouble() != 1.0 or not _replacement_junction_ok(mol, specs, joins[0]):
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
    locants = _assembly_numbering(graph, rings, joins[0], marked, entries, specs, cite_marked, mol, stereo)
    replacement = _assembly_replacement(mol, specs, locants)
    lead = _assembly_stereo_prefix(locants, stereo) if stereo else ""

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
        if any(order not in SUFFIX_OF_ORDER for order in orders.values()):
            return None
        pieces = []
        for order in sorted(set(orders.values())):
            group = [a for a in marked if orders[a] == order]
            spots = ",".join(cite(c) for c in sorted((locants[m] for m in group), key=_locant_order))
            pieces.append(f"{spots}-{valence_word(len(group), order)}")
        core = f"{replacement}[{base}]-" + "-".join(pieces)
        return lead + assembly_join(prefix, core)
    spots = ",".join(cite(c) for c in sorted((locants[m] for m in marked), key=_locant_order))
    core = f"{replacement}[{base}]-{spots}-{multiplied_word(len(frees), 'yl')}"
    return lead + assembly_join(prefix, core)
