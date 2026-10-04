"""Multiplicative names of acyclic structures whose identical chain parents are joined by a heteroatom or a
concatenated linker (P-15.3.1.2, P-15.3.2, P-29.4.2, P-29.5.2): 2,2'-[ethane-1,2-diylbis(oxy)]di(ethan-1-ol),
2,2'-[oxybis(ethane-2,1-diyloxy)]diacetic acid, 3,3'-(methylazanediyl)dipropanoic acid."""

from itertools import combinations

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._hetero_prefixes import is_functional_carbon
from ._multiplicative import _Context
from ._multiplicative_linker import DecompositionRejected, name_component
from ._multiplicative_text import enclose, multiplier_word, unit_phrase
from ._substituents import name_branch

_LINKER_ELEMENTS = {5, 7, 8, 14, 15, 16, 32, 33, 34, 52}
_RUN_ELEMENTS = {8, 16, 34, 52, 14}
_MAX_UNITS = 6
_SKELETAL_UNITS = 4


def _arm(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _key(mol, atoms, attach):
    from ._polyfunctional import _unit_molecule

    unit, _ = _unit_molecule(mol, atoms, attach)
    return Chem.MolToSmiles(unit, isomericSmiles=False)


def _is_linker_atom(mol, atom):
    if atom.GetAtomicNum() not in _LINKER_ELEMENTS or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope():
        return False
    if atom.GetDegree() < 2 or atom.GetNumRadicalElectrons():
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondTypeAsDouble() != 1.0:
            return False
        if other.GetAtomicNum() == 6 and is_functional_carbon(mol, other.GetIdx()):
            return False
    return True


def _principal_atoms(mol):
    from ._polyfunctional import _SENIORITY, _group_of, _ring_occurrences

    groups = {}
    for atom in mol.GetAtoms():
        found = _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], {})[atom.GetIdx()] = found[1]
    ring_groups = _ring_occurrences(mol)
    classes = set(groups) | {c for c, _, _ in ring_groups}
    principal = next((name for name in _SENIORITY if name in classes), None)
    if principal is None:
        return None
    atoms = set(groups.get(principal, {}))
    for owned in groups.get(principal, {}).values():
        atoms |= owned
    for cls, _, owned in ring_groups:
        if cls == principal:
            atoms |= set(owned)
    return principal, atoms


def _ring_system_of(mol, atom):
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    system = next((r for r in rings if atom in r), None)
    if system is None:
        return None
    system = set(system)
    grown = True
    while grown:
        grown = False
        for ring in rings:
            if ring & system and not ring <= system:
                system |= ring
                grown = True
    return system


def _span(mol, graph, handles):
    start = handles[0]
    parent, stack = {start: None}, [start]
    while stack:
        a = stack.pop()
        for n in graph[a]:
            if n not in parent:
                parent[n] = a
                stack.append(n)
    span = set(handles)
    for h in handles[1:]:
        node = h
        while node is not None:
            span.add(node)
            node = parent[node]
    for atom in list(span):
        system = _ring_system_of(mol, atom)
        if system:
            span |= system
    return span


def _components(mol, graph, span):
    comp_of, components = {}, {}
    for atom in sorted(span):
        if atom in comp_of:
            continue
        z = mol.GetAtomWithIdx(atom).GetAtomicNum()
        system = _ring_system_of(mol, atom)
        if system:
            members = system
        else:
            members, stack = {atom}, [atom]
            while stack:
                a = stack.pop()
                for n in graph[a]:
                    neighbor = mol.GetAtomWithIdx(n)
                    if (
                        n in span
                        and n not in members
                        and not neighbor.IsInRing()
                        and neighbor.GetAtomicNum() == z
                        and (z == 6 or z in _RUN_ELEMENTS)
                    ):
                        members.add(n)
                        stack.append(n)
        cid = len(components)
        kind = "ring" if system else "carbon" if z == 6 else "hetero"
        components[cid] = (kind, sorted(members))
        for a in members:
            comp_of[a] = cid
    return _merge_assemblies(mol, graph, comp_of, components)


def _merge_assemblies(mol, graph, comp_of, components):
    parent = {cid: cid for cid in components}

    def find(c):
        while parent[c] != c:
            c = parent[c]
        return c

    for atom, cid in comp_of.items():
        for n in graph[atom]:
            other = comp_of.get(n)
            if other is not None and other != cid and components[cid][0] == "ring" and components[other][0] == "ring":
                parent[find(cid)] = find(other)
    groups = {}
    for cid in components:
        groups.setdefault(find(cid), []).append(cid)
    merged, new_of = {}, {}
    for new, (_, members) in enumerate(sorted(groups.items())):
        atoms = sorted(a for c in members for a in components[c][1])
        kind = "assembly" if len(members) > 1 else components[members[0]][0]
        merged[new] = (kind, atoms)
        for a in atoms:
            new_of[a] = new
    return new_of, merged


def _longest_hetero_run(graph, comp_of, components):
    chain_comps = {cid for cid, (kind, _) in components.items() if kind in ("carbon", "hetero")}
    seen, most = set(), 0
    for start in chain_comps:
        if start in seen:
            continue
        group, stack = {start}, [start]
        while stack:
            cid = stack.pop()
            for atom in components[cid][1]:
                for n in graph[atom]:
                    other = comp_of.get(n)
                    if other in chain_comps and other not in group:
                        group.add(other)
                        stack.append(other)
        seen |= group
        most = max(most, sum(1 for cid in group if components[cid][0] == "hetero"))
    return most


def _is_fused(mol, kind, atoms):
    return kind == "ring" and sum(1 for r in mol.GetRingInfo().AtomRings() if set(r) <= set(atoms)) > 1


def _fused_central(mol, graph, atoms, cedges, arm_atoms):
    from ._diester_ring_diyl import evaluate_skeleton
    from ._multiplicative_linker import Part
    from ._polyfunctional import _require_mancude_system

    ring_set = set(atoms)
    _require_mancude_system(mol, ring_set)
    members = [list(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= ring_set]
    blocked = set().union(*(arm_atoms(graph, y, x) for x, y, _ in cedges))
    try:
        found = evaluate_skeleton(mol, graph, "ring", members, ring_set, [x for x, _, _ in cedges], blocked, "yl")
    except UnsupportedStructure:
        return None
    return Part(found[1], False, True) if found else None


def _component_part(mol, graph, kind, atoms, attachments, ctx, directed):
    if kind != "assembly":
        return name_component(mol, kind, atoms, attachments, ctx, directed)
    if directed is not None or any(order != 1 for _, _, order in attachments):
        raise DecompositionRejected("a ring assembly inside a concatenated linker is not supported")
    from ._chain_assembly import assembly_diyl
    from ._multiplicative_linker import Part

    name = assembly_diyl(mol, graph, halogen_substituents(mol), frozenset(), atoms, [(x, y) for x, y, _ in attachments])
    if name is None:
        raise DecompositionRejected("this ring assembly is not supported as a linker")
    return Part(name, False, True)


def _make_context(mol, graph):
    halogens = halogen_substituents(mol)
    aromatic = frozenset()

    def entry(_mol, owner, root, _ctx):
        return name_branch(graph, root, owner, halogens, aromatic, mol=mol, unsaturated=True)

    return _Context([], None, None, [], set(), {}, entry)


def _linker_text(count, central, arm_parts):
    if not arm_parts:
        return enclose(central.text) if (central.has_prefix or central.has_locants) else central.text
    central_text = enclose(central.text) if central.has_prefix else central.text
    pieces = [enclose(p.text) if p.has_prefix and len(arm_parts) > 1 else p.text for p in arm_parts]
    arm_text = "".join(pieces)
    plain = len(arm_parts) == 1 and arm_parts[0].has_locants and not arm_parts[0].has_prefix
    simple = (plain and not arm_text.startswith(("di", "tri", "tetra"))) or arm_text == "nitrilo"
    arm_enclosed = arm_text if arm_text == "nitrilo" else enclose(arm_text)
    return enclose(central_text + multiplier_word(count, use_bis=not simple) + arm_enclosed)


def _attempt(mol, graph, stereo, arms):
    from ._polyfunctional import _arm_atoms, _select, _stereo_arms, _unit_molecule

    count = len(arms)
    roots = [r for r, _, _ in arms]
    handles = [h for _, h, _ in arms]
    arm_atoms = set().union(*(a for _, _, a in arms))
    span = _span(mol, graph, handles)
    if span & arm_atoms:
        return None
    comp_of, components = _components(mol, graph, span)
    if _longest_hetero_run(graph, comp_of, components) >= _SKELETAL_UNITS:
        return None
    edges = {cid: [] for cid in components}
    for a in span:
        for n in graph[a]:
            order = int(mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble())
            if n in span and comp_of[n] != comp_of[a]:
                edges[comp_of[a]].append((a, n, order))
            elif n in roots:
                edges[comp_of[a]].append((a, n, order))
    center = None
    for cid, cedges in edges.items():
        branches = len(cedges)
        if branches < 2 or count % branches:
            continue
        per_branch = count // branches
        keys = set()
        for x, y, order in cedges:
            branch = _arm_atoms(graph, y, x)
            if sum(1 for r in roots if r in branch) != per_branch:
                keys = None
                break
            keys.add((order, _key(mol, branch, y)))
        if keys is not None and len(keys) == 1:
            center = cid
            break
    if center is None:
        return None
    if stereo:
        for _, where, _ in stereo:
            touched = {where} if isinstance(where, int) else set(where)
            if not touched <= arm_atoms:
                return None
    ctx = _make_context(mol, graph)
    kind, atoms = components[center]
    branches = len(edges[center])
    per_branch = count // branches
    chain = []
    x, y, order = edges[center][0]
    previous, current = x, y
    while current not in roots:
        cid = comp_of[current]
        onward = [e for e in edges[cid] if not (e[0] == current and e[1] == previous)]
        fan = len(onward) == per_branch and all(e[1] in roots for e in onward)
        if not fan and len(onward) != 1:
            return None
        chain.append((cid, (current, previous, order), onward, fan))
        if fan:
            break
        previous, current, order = onward[0][0], onward[0][1], onward[0][2]
    try:
        central = _fused_central(mol, graph, atoms, edges[center], _arm_atoms) if _is_fused(mol, kind, atoms) else None
        if central is None:
            central = _component_part(mol, graph, kind, atoms, edges[center], ctx, None)
        arm_parts = []
        for cid, toward_center, onward, fan in chain:
            akind, aatoms = components[cid]
            if fan and per_branch > 1:
                arm_parts.append(_component_part(mol, graph, akind, aatoms, [toward_center, *onward], ctx, None))
            else:
                arm_parts.append(
                    _component_part(
                        mol, graph, akind, aatoms, [toward_center, onward[0]], ctx, (onward[0][0], toward_center[0])
                    )
                )
    except DecompositionRejected:
        return None
    parts_atoms = [a for _, _, a in arms]
    units = [_unit_molecule(mol, atoms_, root) for atoms_, root in zip(parts_atoms, roots)]
    stereo_text = ""
    if stereo:
        parts, stereo_text = _stereo_arms(mol, stereo, parts_atoms, units)
    else:
        unit, attach = units[0]
        try:
            _, _, parts = _select(unit, attach)
        except UnsupportedStructure:
            return None
    prefix, body, tail, locant = parts[:4]
    lead = ",".join(str(locant) + "'" * i for i in range(count)) + "-" if locant is not None else ""
    linker = _linker_text(branches, central, arm_parts)
    text = prefix + body
    if prefix:
        if tail.startswith(" "):
            return f"{stereo_text}{lead}{linker}{multiplier_word(count, True)}({text}{tail})"
        return f"{stereo_text}{lead}{linker}{multiplier_word(count, True)}({text}){tail}"
    return f"{stereo_text}{lead}{linker}{multiplier_word(count, False)}{unit_phrase(text, tail)}"


def chain_multiplicative_name(mol, stereo):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    graph = adjacency(mol)
    found = _principal_atoms(mol)
    if found is None:
        return None
    principal, anchors = found
    linkers = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if _is_linker_atom(mol, a) and not (principal == "amine" and a.GetAtomicNum() == 7)
    }
    candidates = {}
    for r in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(r).GetAtomicNum() != 6 or mol.GetAtomWithIdx(r).IsInRing():
            continue
        for h in graph[r]:
            if h in linkers or (mol.GetAtomWithIdx(h).IsInRing() and mol.GetBondBetweenAtoms(r, h).GetBondTypeAsDouble() == 1.0):
                atoms = _arm(graph, r, h)
                if h not in atoms:
                    candidates.setdefault(_key(mol, atoms, r), []).append((r, h, atoms))
    ranked = []
    for key, arms in candidates.items():
        if len(arms) < 2:
            continue
        for size in range(min(len(arms), _MAX_UNITS), 1, -1):
            for subset in combinations(arms, size):
                if any(a[2] & b[2] for a, b in combinations(subset, 2)):
                    continue
                if not anchors <= set().union(*(a[2] for a in subset)) or any(not anchors & a[2] for a in subset):
                    continue
                ranked.append((-size, key, list(subset)))
    ranked.sort(key=lambda item: item[:2])
    for _, _, arms in ranked:
        name = _attempt(mol, graph, stereo, arms)
        if name is not None:
            return name
    return None


def has_chain_multiplicative_shape(mol):
    try:
        return chain_multiplicative_name(mol, []) is not None
    except Exception:
        return False
