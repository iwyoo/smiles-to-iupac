"""Multiplicative names of acyclic structures whose identical chain parents are joined by a heteroatom or a
concatenated linker (P-15.3.1.2, P-15.3.2, P-29.4.2, P-29.5.2): 2,2'-[ethane-1,2-diylbis(oxy)]di(ethan-1-ol),
2,2'-[oxybis(ethane-2,1-diyloxy)]diacetic acid, 3,3'-(methylazanediyl)dipropanoic acid."""

import re
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
_SUBSTITUTED_PREFIX = re.compile(r"(?:carboxy|hydroxy|amino|chloro|bromo|fluoro|iodo|cyano|oxo|nitro|sulfanyl|methoxy|ethoxy)[a-z]+")
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
    heads = set(atoms)
    for owned in groups.get(principal, {}).values():
        atoms |= owned
    for cls, _, owned in ring_groups:
        if cls == principal:
            atoms |= set(owned)
            heads.add(min(owned))
    return principal, atoms, heads


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
    return _grow_assembly(mol, graph, span)


def _grow_assembly(mol, graph, span):
    """Identical ring systems joined directly to the linker's rings belong to the same ring assembly (P-28.2)."""
    from ._multiplicative import _bare_key

    grown = True
    while grown:
        grown = False
        for atom in sorted(span):
            system = _ring_system_of(mol, atom)
            if not system:
                continue
            for member in system:
                for n in graph[member]:
                    if n in span or not mol.GetAtomWithIdx(n).IsInRing() or mol.GetBondBetweenAtoms(member, n).GetBondTypeAsDouble() != 1.0:
                        continue
                    other = _ring_system_of(mol, n)
                    if other and not other & span and _bare_key(mol, other) == _bare_key(mol, system):
                        span |= other
                        grown = True
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


def _ylidene_ring_part(mol, graph, atoms, attachments):
    """naphthalene-2,3-diylidene: the ketone name of the ring with each ylidene bond replaced by C=O (P-29.3.4.2)."""
    from ._multiplicative_linker import Part
    from ._multiplicative_prefix import probe_name

    external = {y for _, y, _ in attachments}
    keep, stack = set(atoms), list(atoms)
    while stack:
        for n in graph[stack.pop()]:
            if n not in keep and n not in external:
                keep.add(n)
                stack.append(n)
    editable = Chem.RWMol(mol)
    for x, y, _ in attachments:
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(x, oxygen, Chem.BondType.DOUBLE)
        keep.add(oxygen)
    for idx in sorted(set(range(editable.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(idx)
    probe = editable.GetMol()
    Chem.SanitizeMol(probe)
    name = probe_name(Chem.MolToSmiles(probe))
    match = re.fullmatch(r"(.+)-(di|tri|tetra)one", name)
    if match is None:
        raise DecompositionRejected(f"the ylidene ring {name!r} is not a ketone name")
    head, count = match.groups()
    return Part(f"{head}-{count}ylidene", head[0].isdigit() or head[0] in "([", True)


def _component_part(mol, graph, kind, atoms, attachments, ctx, directed):
    if kind == "ring" and attachments and all(order == 2 for _, _, order in attachments):
        if directed is not None:
            raise DecompositionRejected("a ylidene ring inside a concatenated linker is not supported")
        return _ylidene_ring_part(mol, graph, atoms, attachments)
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
        name, compound = name_branch(graph, root, owner, halogens, aromatic, mol=mol, unsaturated=True)
        if (compound or _SUBSTITUTED_PREFIX.search(name)) and not any(ch.isdigit() for ch in name):
            return f"({name})", False
        return name, compound

    return _Context([], None, None, [], set(), {}, entry)


def _linker_text(count, central, arm_parts, led=True):
    if not arm_parts:
        return enclose(central.text) if (central.has_prefix or central.has_locants) and led else central.text
    central_text = enclose(central.text) if central.has_prefix else central.text
    pieces = [enclose(p.text) if p.has_prefix and len(arm_parts) > 1 else p.text for p in arm_parts]
    arm_text = "".join(pieces)
    plain = len(arm_parts) == 1 and arm_parts[0].has_locants and not arm_parts[0].has_prefix
    simple = (plain and not arm_text.startswith(("di", "tri", "tetra"))) or arm_text == "nitrilo"
    arm_enclosed = arm_text if arm_text == "nitrilo" else enclose(arm_text)
    return enclose(central_text + multiplier_word(count, use_bis=not simple) + arm_enclosed)


def _attempt(mol, graph, stereo, arms, unit_kind="chain"):
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
        ylidene = all(order == 2 for _, _, order in edges[center])
        central = (
            _fused_central(mol, graph, atoms, edges[center], _arm_atoms)
            if _is_fused(mol, kind, atoms) and not ylidene
            else None
        )
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
    if unit_kind != "chain":
        if stereo:
            return None
        text = _unit_name_by_pipeline(units[0][0], unit_kind)
        if text is None:
            return None
        lead = ""
        if unit_kind == "amide":
            lead = ",".join("N" + "'" * i for i in range(count)) + "-"
        elif _MULTIPLIED_HYDRIDE.match(text):
            lead = ",".join("1" + "'" * i for i in range(count)) + "-"
        linker = _linker_text(branches, central, arm_parts, bool(lead) or unit_kind != "hydride")
        if unit_kind == "hydride" or text.startswith(("N-", "di", "tri", "tetra")) or any(ch.isdigit() or ch == "-" for ch in text):
            return f"{lead}{linker}{multiplier_word(count, True)}({text})"
        return f"{lead}{linker}{multiplier_word(count, False)}{text}"
    linker = _linker_text(branches, central, arm_parts)
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
    text = prefix + body
    if prefix:
        if tail.startswith(" "):
            return f"{stereo_text}{lead}{linker}{multiplier_word(count, True)}({text}{tail})"
        return f"{stereo_text}{lead}{linker}{multiplier_word(count, True)}({text}){tail}"
    return f"{stereo_text}{lead}{linker}{multiplier_word(count, False)}{unit_phrase(text, tail)}"


_HYDRIDE_ENDINGS = ("phosphane", "arsane", "silane", "germane", "stannane", "plumbane")
_HYDRIDE_ORDER = (15, 33, 14, 32, 50, 82)
_MULTIPLIED_HYDRIDE = re.compile(r"^(?:di|tri|tetra|penta|hexa|hepta|octa)(?:phosphane|arsane|silane|germane|stannane|plumbane)$")


def _unit_name_by_pipeline(unit, unit_kind):
    from .core import smiles_to_iupac

    editable = Chem.RWMol(unit)
    for atom in editable.GetAtoms():
        if atom.GetAtomMapNum() and atom.GetNoImplicit():
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(editable)
    try:
        name = smiles_to_iupac(Chem.MolToSmiles(editable))
    except UnsupportedStructure:
        return None
    if unit_kind == "amide":
        return name if name.endswith("amide") else None
    return name if name.endswith(_HYDRIDE_ENDINGS) else None


def _rank_candidates(candidates, anchors, heads=None, mol=None):
    ranked = []
    for key, arms in candidates.items():
        if len(arms) < 2:
            continue
        for size in range(min(len(arms), _MAX_UNITS), 1, -1):
            for subset in combinations(arms, size):
                if any(a[2] & b[2] for a, b in combinations(subset, 2)):
                    continue
                if anchors is not None:
                    covered = set().union(*(a[2] for a in subset))
                    if any(not anchors & a[2] for a in subset):
                        continue
                    if not anchors <= covered and not _beats_substitutive(mol, heads, covered):
                        continue
                ranked.append((-size, key, list(subset)))
    ranked.sort(key=lambda item: item[:2])
    return [arms for _, _, arms in ranked]


def _beats_substitutive(mol, heads, covered):
    from ._polyfunctional import _select

    if mol is None or heads is None:
        return False
    try:
        key, _, _ = _select(mol)
    except UnsupportedStructure:
        return True
    return isinstance(key[0], int) and len(heads & covered) > -key[0]


def _bare_chain(graph, root, senior_atoms, arm_atoms):
    if arm_atoms != senior_atoms:
        return False
    end = [a for a in senior_atoms if sum(n in senior_atoms for n in graph[a]) <= 1]
    return root in end and all(sum(n in senior_atoms for n in graph[a]) <= 2 for a in senior_atoms)


def _hydride_candidates(mol, graph):
    if any(a.GetAtomicNum() == 7 for a in mol.GetAtoms()):
        return {}
    present = {a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() in _HYDRIDE_ORDER}
    element = next((z for z in _HYDRIDE_ORDER if z in present), None)
    if element is None:
        return {}
    senior = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == element}
    candidates = {}
    for r in senior:
        atom = mol.GetAtomWithIdx(r)
        if atom.IsInRing() or atom.GetFormalCharge():
            continue
        for h in graph[r]:
            if h in senior or mol.GetAtomWithIdx(h).GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(r, h).GetBondTypeAsDouble() != 1.0:
                continue
            atoms = _arm(graph, r, h)
            inside = senior & atoms
            if h in atoms or (inside != {r} and not _bare_chain(graph, r, inside, atoms)):
                continue
            candidates.setdefault(_key(mol, atoms, r), []).append((r, h, atoms))
    return candidates


def chain_multiplicative_name(mol, stereo):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    graph = adjacency(mol)
    found = _principal_atoms(mol)
    if found is None:
        hydrides = _hydride_candidates(mol, graph)
        if not hydrides:
            return None
        candidates = hydrides
        for arms in _rank_candidates(candidates, None):
            name = _attempt(mol, graph, stereo, arms, "hydride")
            if name is not None:
                return name
        return None
    principal, anchors, heads = found
    linkers = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if _is_linker_atom(mol, a) and not (principal == "amine" and a.GetAtomicNum() == 7)
    }
    candidates, nitrogen = {}, {}
    for r in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(r)
        if atom.IsInRing():
            continue
        for h in graph[r]:
            bond = mol.GetBondBetweenAtoms(r, h)
            if atom.GetAtomicNum() == 6:
                ok = h in linkers or (mol.GetAtomWithIdx(h).IsInRing() and bond.GetBondTypeAsDouble() in (1.0, 2.0))
                target = candidates
            elif atom.GetAtomicNum() == 7 and principal in ("amide", "sulfonamide") and r in anchors:
                other = mol.GetAtomWithIdx(h)
                ok = bond.GetBondTypeAsDouble() == 1.0 and other.GetAtomicNum() == 6 and not is_functional_carbon(mol, h)
                target = nitrogen
            else:
                continue
            if ok:
                atoms = _arm(graph, r, h)
                if h not in atoms:
                    target.setdefault(_key(mol, atoms, r), []).append((r, h, atoms))
    for arms in _rank_candidates(candidates, anchors, heads, mol):
        name = _attempt(mol, graph, stereo, arms)
        if name is not None:
            return name
    for arms in _rank_candidates(nitrogen, anchors, heads, mol):
        name = _attempt(mol, graph, stereo, arms, "amide")
        if name is not None:
            return name
    return None


def has_chain_multiplicative_shape(mol):
    try:
        return chain_multiplicative_name(mol, []) is not None
    except Exception:
        return False
