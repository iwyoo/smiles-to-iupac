"""Multiplicative names of acyclic structures whose identical chain parents are joined by a heteroatom or a
concatenated linker (P-15.3.1.2, P-15.3.2, P-29.4.2, P-29.5.2): 2,2'-[ethane-1,2-diylbis(oxy)]di(ethan-1-ol),
2,2'-[oxybis(ethane-2,1-diyloxy)]diacetic acid, 3,3'-(methylazanediyl)dipropanoic acid."""

import re
from itertools import combinations

from rdkit import Chem

from ._chalcogen_chain_heterone import chalcogen_chain_oxygen_count
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents, is_prefix_placeholder
from ._diacyl_linker import acyl_part, diacyl_linker, diamido_linker
from ._hetero_prefixes import is_functional_carbon
from ._multiplicative import _Context
from ._multiplicative_groups import terminal_ligand
from ._multiplicative_linker import DecompositionRejected, name_component
from ._multiplicative_text import enclose, join_parts, multiplier_word, opens_with_assembly, unit_phrase
from ._substituents import name_branch

_LINKER_ELEMENTS = {5, 7, 8, 14, 15, 16, 32, 33, 34, 50, 51, 52, 82, 83}
_RUN_ELEMENTS = {7, 8, 14, 15, 16, 32, 33, 34, 50, 51, 52, 82, 83}
_MAX_UNITS = 6
_SUBSTITUTED_PREFIX = re.compile(r"(?:carboxy|hydroxy|amino|chloro|bromo|fluoro|iodo|cyano|oxo|nitro(?!so\b)|sulfanyl|methoxy|ethoxy)[a-z]+")
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

    try:
        unit, _ = _unit_molecule(mol, atoms, attach)
    except UnsupportedStructure:
        return ("arm", tuple(sorted(atoms)), attach)
    return Chem.MolToSmiles(unit, isomericSmiles=False)


_DOUBLY_BONDED_PAIRS = (7, 14, 15, 32, 33, 50, 51, 82, 83)


def _acyl_carbon(carbon):
    return any(
        b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(carbon).GetAtomicNum() == 8 and b.GetOtherAtom(carbon).GetDegree() == 1
        for b in carbon.GetBonds()
    )


def _acyl_branch(graph, carbon, nitrogen):
    """The atoms of the acyl group on `carbon` seen from `nitrogen`, or none when they close a ring back to it."""
    atoms = _arm(graph, carbon, nitrogen)
    return set() if nitrogen in atoms else atoms


def _is_linker_atom(mol, atom, acyl_prefixes=frozenset()):
    """`acyl_prefixes`: the carbons that a nitrogen linker may carry as acyl groups, cited as prefixes of the linker
    because a group senior to the amide is the principal one (P-41, P-65.1)."""
    if atom.GetAtomicNum() not in _LINKER_ELEMENTS or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope():
        return False
    if atom.GetDegree() < 2 or atom.GetNumRadicalElectrons():
        return False
    imine_nitrogen = (
        atom.GetAtomicNum() == 7
        and atom.GetDegree() == 2
        and sorted(b.GetBondTypeAsDouble() for b in atom.GetBonds()) == [1.0, 2.0]
        and all(b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds())
    )
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        azo = atom.GetAtomicNum() == other.GetAtomicNum() and atom.GetAtomicNum() in _DOUBLY_BONDED_PAIRS and bond.GetBondTypeAsDouble() == 2.0
        sulfonyl_oxygen = (
            atom.GetAtomicNum() == 16
            and bond.GetBondTypeAsDouble() == 2.0
            and (terminal_ligand(other, atom) or (None,))[0] in ("O", "NH")
        )
        if bond.GetBondTypeAsDouble() != 1.0 and not azo and not imine_nitrogen and not sulfonyl_oxygen:
            return False
        if other.GetAtomicNum() == 6 and is_functional_carbon(mol, other.GetIdx()) and not imine_nitrogen:
            acyl_prefix = atom.GetAtomicNum() == 7 and other.GetIdx() in acyl_prefixes
            if not acyl_prefix and not _diacyl_carbon(mol, other, atom):
                return False
    return True


def _diacyl_carbon(mol, carbon, linker):
    """A carbonyl carbon bonded to the `linker` heteroatom whose carbon skeleton, side groups allowed, reaches
    exactly one second carbonyl carbon bonded to a heteroatom (a diacyl linking group)."""
    hetero = (7, 8, 16, 34, 52)

    def terminals(atom):
        return [
            n.GetIdx()
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() in hetero
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]

    if linker.GetAtomicNum() not in hetero or len(terminals(carbon)) != 1:
        return False
    seen, stack, ends = {linker.GetIdx(), carbon.GetIdx()}, [carbon], 0
    while stack:
        current = stack.pop()
        for n in current.GetNeighbors():
            if n.GetIdx() in seen or n.GetIdx() in terminals(current):
                continue
            seen.add(n.GetIdx())
            if n.GetAtomicNum() != 6:
                if current is carbon:
                    return n.GetAtomicNum() in hetero and ends == 0 and not stack
                continue
            if terminals(n):
                outward = [m for m in n.GetNeighbors() if m.GetIdx() not in seen and m.GetIdx() not in terminals(n)]
                if len(terminals(n)) != 1 or len(outward) != 1 or outward[0].GetAtomicNum() not in hetero:
                    return False
                ends += 1
                continue
            stack.append(n)
    return ends == 1


def _principal_atoms(mol):
    from ._polyfunctional import _SENIORITY, _group_of, _ring_occurrences, _substituted_amine_nitrogen

    groups = {}
    for atom in mol.GetAtoms():
        found = _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], {})[atom.GetIdx()] = found[1]
        elif _substituted_amine_nitrogen(mol, atom) and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors()):
            groups.setdefault("amine", {})[atom.GetIdx()] = {atom.GetIdx()}
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
    from ._appendix3_skeletons import appendix3_multivalent_group
    from ._polyfunctional import _require_mancude_system

    ring_set = set(atoms)
    natural_product = appendix3_multivalent_group(mol, graph, ring_set, [(x, y) for x, y, _ in cedges])
    if natural_product is not None:
        return Part(natural_product[0], natural_product[1], True)
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


def _mixed_ring_part(mol, graph, atoms, attachments):
    from ._diester_ring_diyl import evaluate_skeleton
    from ._multiplicative_linker import Part
    from ._polyfunctional import _arm_atoms, _require_mancude_system

    ring_set = set(atoms)
    members = [list(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= ring_set]
    if len(members) > 1:
        _require_mancude_system(mol, ring_set)
    blocked = set().union(*(_arm_atoms(graph, y, x) for x, y, _ in attachments))
    orders = {x: order for x, _, order in attachments}
    try:
        found = evaluate_skeleton(mol, graph, "ring", members, ring_set, list(orders), blocked, "yl", orders=orders)
    except UnsupportedStructure:
        return None
    return Part(found[1], False, True) if found else None


def _component_part(mol, graph, kind, atoms, attachments, ctx, directed, first=None, stereo=None):
    mixed = len({order for _, _, order in attachments}) > 1
    if kind == "ring" and mixed:
        if directed is not None:
            raise DecompositionRejected("a yl-ylidene ring inside a concatenated linker is not supported")
        part = _mixed_ring_part(mol, graph, atoms, attachments)
        if part is None:
            raise DecompositionRejected("this yl-ylidene ring is not supported")
        return part
    if kind == "ring" and attachments and all(order == 2 for _, _, order in attachments):
        if directed is not None:
            raise DecompositionRejected("a ylidene ring inside a concatenated linker is not supported")
        return _ylidene_ring_part(mol, graph, atoms, attachments)
    if kind != "assembly":
        return name_component(mol, kind, atoms, attachments, ctx, directed, first, stereo or ())
    if directed is not None or any(order != 1 for _, _, order in attachments) and not mixed:
        raise DecompositionRejected("a ring assembly inside a concatenated linker is not supported")
    from ._chain_assembly import assembly_diyl
    from ._multiplicative_linker import Part

    name = assembly_diyl(mol, graph, halogen_substituents(mol), frozenset(), atoms, [(x, y) for x, y, _ in attachments], stereo)
    if name is None:
        raise DecompositionRejected("this ring assembly is not supported as a linker")
    return Part(name, False, True)


def _context_element(mol, element):
    kind, where, code = element
    if kind == "bond":
        return kind, mol.GetBondBetweenAtoms(*where).GetIdx(), code
    return element


def _make_context(mol, graph, stereo=()):
    halogens = halogen_substituents(mol)
    aromatic = frozenset()

    def entry(_mol, owner, root, _ctx):
        name, compound = name_branch(graph, root, owner, halogens, aromatic, mol=mol, unsaturated=True)
        if (compound or _SUBSTITUTED_PREFIX.search(name)) and not any(ch.isdigit() for ch in name):
            return f"({name})", False
        return name, compound

    return _Context([], None, None, list(stereo), set(), {}, entry)


_MULTIPLIED_HETERO_HYDRIDE = re.compile(
    r"(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)(?:silane|germane|stannane|plumbane|phosphane|arsane|stibane|bismuthane|azane|sulfane|selane|tellane|borane|alumane)"
)


def _linker_text(count, central, arm_parts, led=True, bare=False):
    if not arm_parts:
        enclosed = (central.has_prefix or central.has_locants) and led and not (bare and opens_with_assembly(central.text))
        return enclose(central.text) if enclosed else central.text
    central_text = enclose(central.text) if central.has_prefix else central.text
    pieces = [enclose(p.text) if p.has_prefix and len(arm_parts) > 1 else p.text for p in arm_parts]
    arm_text = join_parts(pieces)
    plain = len(arm_parts) == 1 and arm_parts[0].has_locants and not arm_parts[0].has_prefix
    simple = (plain and not _MULTIPLIED_HETERO_HYDRIDE.match(arm_text)) or arm_text == "nitrilo"
    arm_enclosed = arm_text if arm_text == "nitrilo" else enclose(arm_text)
    return enclose(central_text + multiplier_word(count, use_bis=not simple) + arm_enclosed)


def _skeletal_linker_name(mol, graph, arms, span):
    """Two units joined by a chain that is itself a skeletal replacement parent: the linker is cited by its 'a' name
    (P-51.3.1)."""
    from ._heteroacyclic import skeletal_linker
    from ._polyfunctional import _select, _unit_molecule

    roots = [r for r, _, _ in arms]
    handles = [h for _, h, _ in arms]
    linker = skeletal_linker(mol, graph, span, handles, frozenset(roots))
    if linker is None:
        return None
    unit, attach = _unit_molecule(mol, arms[0][2], roots[0])
    try:
        _, _, parts = _select(unit, attach)
    except UnsupportedStructure:
        return None
    prefix, body, tail, locant = parts[:4]
    lead = ",".join(str(locant) + "'" * i for i in range(2)) + "-" if locant is not None else ""
    text = prefix + body
    linker = enclose(linker)
    if prefix:
        return f"{lead}{linker}{multiplier_word(2, True)}({text}){tail}"
    return f"{lead}{linker}{multiplier_word(2, False)}{unit_phrase(text, tail)}"


def _bare_groups_on_one_ring_atom(mol, kind, edges, arms):
    """Principal groups with one carbon each on a ring are a multiplied suffix of the ring parent, not units of a
    multiplicative name (P-65.1.2.2.2, P-66.1.1.1.1.3)."""
    return (
        kind == "ring"
        and all(sum(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in atoms) == 1 for _, _, atoms in arms)
    )


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
    # P-21.2.3.1: where nitrogen is present the amine name outranks skeletal replacement
    if unit_kind != "amine" and _longest_hetero_run(graph, comp_of, components) >= _SKELETAL_UNITS:
        if unit_kind != "chain" or count != 2 or stereo:
            return None
        return _skeletal_linker_name(mol, graph, arms, span)
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
            keys.add((_key(mol, branch, y), order))
        if keys is not None and len(keys) == 1:
            center = cid
            break
    if center is None:
        return None
    unit_stereo, branch_stereo = [], [[] for _ in edges[center]]
    branch_atoms = [_arm_atoms(graph, y, x) for x, y, _ in edges[center]]
    kind, atoms = components[center]
    centre_stereo = {}
    ring_stereo = []
    for element in stereo:
        touched = {element[1]} if isinstance(element[1], int) else set(element[1])
        if touched <= arm_atoms:
            unit_stereo.append(element)
            continue
        if kind == "assembly" and element[0] == "atom" and touched <= set(atoms):
            centre_stereo[element[1]] = element[2]
            continue
        if kind == "ring" and element[0] == "atom" and touched <= set(atoms):
            ring_stereo.append(element)
            continue
        holders = [i for i, atoms_ in enumerate(branch_atoms) if touched <= atoms_]
        if len(holders) != 1:
            return None
        branch_stereo[holders[0]].append(element)
    if len({tuple(sorted(code for _, _, code in found)) for found in branch_stereo}) != 1:
        return None
    ctx = _make_context(mol, graph, [_context_element(mol, e) for e in branch_stereo[0]])
    # P-15.6.1.5, P-51.5: a ring bearing a side chain with a further principal group is a substituent of the chain
    # parent that carries them, not a multiplicative centre
    found = _principal_atoms(mol) if unit_kind == "chain" else None
    if found is not None and not found[2] <= arm_atoms and any(mol.GetAtomWithIdx(a).IsInRing() for a in atoms):
        return None
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
            central = _component_part(
                mol, graph, kind, atoms, edges[center], ctx, None, stereo=centre_stereo or ring_stereo or None
            )
        arm_parts = []
        for cid, toward_center, onward, fan in chain:
            akind, aatoms = components[cid]
            if fan and per_branch > 1:
                arm_parts.append(
                    _component_part(mol, graph, akind, aatoms, [toward_center, *onward], ctx, None, toward_center[0])
                )
            else:
                arm_parts.append(
                    _component_part(
                        mol, graph, akind, aatoms, [toward_center, onward[0]], ctx, (onward[0][0], toward_center[0])
                    )
                )
    except DecompositionRejected:
        return None
    except UnsupportedStructure:
        if len({o for _, _, o in edges[center]}) > 1 or _bare_groups_on_one_ring_atom(mol, kind, edges[center], arms):
            return None
        raise
    if count == 2:
        diacyl = diacyl_linker(
            mol, graph, set(atoms), edges[center], [set(components[comp_of[o]][1]) if o in comp_of else set() for _, o, _ in edges[center]], arm_atoms
        )
        if diacyl is not None:
            central = acyl_part(diacyl[0], diacyl[2])
            arm_parts = arm_parts[diacyl[1]:]
        if unit_kind == "chain":
            diamido = diamido_linker(mol, graph, set(range(mol.GetNumAtoms())) - arm_atoms, handles)
            if diamido is not None:
                central, arm_parts = diamido, []
    parts_atoms = [a for _, _, a in arms]
    units = [_unit_molecule(mol, atoms_, root) for atoms_, root in zip(parts_atoms, roots)]
    if unit_kind != "chain":
        if stereo:
            return None
        chain_unit = _substituted_chain_unit(mol, graph, roots[0], parts_atoms[0]) if unit_kind == "hydride" else None
        if unit_kind == "urea":
            text = _urea_unit_text(mol, graph, roots[0], arms[0][1], parts_atoms[0])
        else:
            text = chain_unit or _unit_name_by_pipeline(units[0][0], unit_kind)
        if text is None:
            return None
        lead = ""
        if unit_kind == "anion":
            linker = _linker_text(branches, central, arm_parts, True)
            return f"{linker}{multiplier_word(count, True)}({text})"
        if unit_kind == "ester":
            from ._hetero_prefixes import _without_enclosing_marks

            linker = _without_enclosing_marks(_linker_text(branches, central, arm_parts))
            compound = not text.isalpha()
            return f"{linker} {multiplier_word(count, compound)}{enclose(text) if compound else text}"
        if unit_kind == "polyamine":
            lead = _polyamine_lead(units[0][0], count)
            if lead is None:
                return None
            return f"{lead}{_linker_text(branches, central, arm_parts)}{multiplier_word(count, False)}{unit_phrase(text, '')}"
        if unit_kind in ("amide", "amine", "imine"):
            lead = ",".join("N" + "'" * i for i in range(count)) + "-"
        elif unit_kind == "urea":
            lead = ",".join("N" + "'" * (2 * i) for i in range(count)) + "-"
        elif chain_unit or _MULTIPLIED_HYDRIDE.match(text):
            lead = ",".join("1" + "'" * i for i in range(count)) + "-"
        linker = _linker_text(branches, central, arm_parts, bare=not lead)
        if unit_kind == "inorganic_amide":
            from ._hetero_prefixes import _without_enclosing_marks

            return f"{_without_enclosing_marks(_linker_text(branches, central, arm_parts, False))}{multiplier_word(count, True)}{enclose(text)}"
        if unit_kind == "imine":
            if _imine_unit_needs_bis(text, units[0][0]):
                return f"{lead}{linker}{multiplier_word(count, True)}({text})"
            return f"{lead}{linker}{multiplier_word(count, False)}{unit_phrase(text, '')}"
        if (
            unit_kind == "hydride"
            or text.startswith(("N-", "N'", "di", "tri", "tetra", "thio", "seleno", "telluro"))
            or any(ch.isdigit() or ch in "- " for ch in text)
        ):
            return f"{lead}{linker}{multiplier_word(count, True)}({text})"
        return f"{lead}{linker}{multiplier_word(count, False)}{text}"
    linker = _linker_text(branches, central, arm_parts)
    stereo_text = ""
    if len(ctx.used) != len(branch_stereo[0]):
        return None
    unit, attach = units[0]
    try:
        _, _, parts = _select(unit, attach)
    except UnsupportedStructure:
        amine = None if stereo or unit_stereo else _amine_arm(unit, attach)
        if amine is None:
            return None
        lead = ",".join(str(amine[1]) + "'" * i for i in range(count)) + "-"
        return f"{lead}{linker}{multiplier_word(count, True)}({amine[0]})"
    if unit_stereo:
        parts, stereo_text = _stereo_arms(mol, unit_stereo, parts_atoms, units)
    prefix, body, tail, locant = parts[:4]
    lead = ",".join(str(locant) + "'" * i for i in range(count)) + "-" if locant is not None else ""
    text = prefix + body
    if prefix:
        if tail.startswith(" "):
            return f"{stereo_text}{lead}{linker}{multiplier_word(count, True)}({text}{tail})"
        return f"{stereo_text}{lead}{linker}{multiplier_word(count, True)}({text}){tail}"
    return f"{stereo_text}{lead}{linker}{multiplier_word(count, False)}{unit_phrase(text, tail)}"


_CENTER_ACIDS = ("phosphonic acid", "arsonic acid", "stibonic acid", "phosphonous acid", "arsonous acid", "stibonous acid")
_HYDRIDE_STEMS = ("phosphane", "arsane", "silane", "germane", "stannane", "plumbane")
_HETERONE_SUFFIXES = ("one", "thione", "selone", "tellone")
_HYDRIDE_ENDINGS = _HYDRIDE_STEMS + tuple(
    (stem[:-1] if suffix[0] in "aeiou" else stem) + suffix for stem in _HYDRIDE_STEMS for suffix in _HETERONE_SUFFIXES
)
_HYDRIDE_ORDER = (15, 33, 14, 32, 50, 82)
_MULTIPLIED_HYDRIDE = re.compile(r"^(?:di|tri|tetra|penta|hexa|hepta|octa)(?:phosphane|arsane|silane|germane|stannane|plumbane)$")


def _imine_unit_needs_bis(text, unit):
    """'bis(1-phenylmethanimine)' when the unit name opens with a locant or a prefix that 'di' would fuse with;
    'di(propan-2-imine)' otherwise (P-16.3.6)."""
    return text[0].isdigit() or text.startswith(("N", "di", "tri", "tetra")) or _ring_on_imine_carbon(unit)


def _ring_on_imine_carbon(unit):
    """A ring joined to the carbon of the C=N group: 'bis(phenylmethanimine)' keeps the prefix from reading as 'diphenyl'."""
    return any(
        a.GetAtomicNum() == 6
        and not a.IsInRing()
        and any(n.IsInRing() for n in a.GetNeighbors())
        and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(a).GetAtomicNum() == 7 for b in a.GetBonds())
        for a in unit.GetAtoms()
    )


def _amine_arm(unit, attach_atom):
    from ._polyfunctional import amine_unit_name

    editable = Chem.RWMol(unit)
    attach = next((a.GetIdx() for a in editable.GetAtoms() if a.GetAtomMapNum() and not is_prefix_placeholder(a)), None)
    if attach is None:
        return None
    for atom in editable.GetAtoms():
        if not is_prefix_placeholder(atom):
            atom.SetAtomMapNum(0)
    mol = editable.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.MolSanitizeException:
        return None
    return amine_unit_name(mol, attach)


def _unit_name_by_pipeline(unit, unit_kind):
    from .core import smiles_to_iupac

    editable = Chem.RWMol(unit)
    for atom in editable.GetAtoms():
        if is_prefix_placeholder(atom):
            continue
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
    if unit_kind == "imine":
        return name if name.endswith("imine") else None
    if unit_kind == "amine":
        single_nitrogen = sum(a.GetAtomicNum() == 7 for a in unit.GetAtoms()) == 1
        return name if single_nitrogen and name.endswith(("amine", "aniline")) else None
    if unit_kind == "polyamine":
        return name if name.endswith("diamine") else None
    if unit_kind == "phosphonic":
        return name if name in _CENTER_ACIDS else None
    if unit_kind == "centre_acid":
        return name if _CENTRE_ACID.search(name) else None
    if unit_kind == "carbamic":
        return name if _CARBAMIC_ACID.search(name) else None
    if unit_kind == "inorganic_amide":
        return name if _INORGANIC_AMIDE.search(name) else None
    if unit_kind == "anion":
        return name if re.search(r"(?:ide|uide|ate)$", name) else None
    if unit_kind == "ester":
        from ._acid_derivatives import anion_name

        return anion_name(name, letters_in_anion=False) if name.endswith(" acid") else None
    return name if name.endswith(_HYDRIDE_ENDINGS) else None


_MARKER = "tert-butyl"
_MARKER_LOCANT = re.compile(r"(N['\u2032]*\d+)-(?:\d+[a-z]*-)?" + _MARKER)


def _polyamine_lead(unit, count):
    """The N locants of the junction nitrogen in the identical parents, `N1,N1'-`: the unit named with a tert-butyl
    group on that nitrogen shows its locant (P-62.2.5.2)."""
    from .core import smiles_to_iupac

    editable = Chem.RWMol(unit)
    attach = next(a.GetIdx() for a in editable.GetAtoms() if a.GetAtomMapNum() and not is_prefix_placeholder(a))
    tail = editable.AddAtom(Chem.Atom(6))
    editable.AddBond(attach, tail, Chem.BondType.SINGLE)
    for _ in range(3):
        methyl = editable.AddAtom(Chem.Atom(6))
        editable.AddBond(tail, methyl, Chem.BondType.SINGLE)
    for atom in editable.GetAtoms():
        if is_prefix_placeholder(atom):
            continue
        if atom.GetAtomMapNum() and atom.GetNoImplicit():
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        atom.SetAtomMapNum(0)
    try:
        Chem.SanitizeMol(editable)
        name = smiles_to_iupac(Chem.MolToSmiles(editable))
    except (UnsupportedStructure, Chem.rdchem.MolSanitizeException):
        return None
    found = _MARKER_LOCANT.search(name)
    if found is None:
        return None
    locant = found.group(1)
    return ",".join(locant + "'" * i for i in range(count)) + "-"


def anion_multiplicative_name(mol, centers):
    """Identical anionic parents joined by a linker, e.g. (1,4-phenylene)bis(phosphanide) (P-72.5.1.1)."""
    if len(Chem.GetMolFrags(mol)) != 1 or len(centers) < 2:
        return None
    graph = adjacency(mol)
    candidates = {}
    for r in centers:
        if mol.GetAtomWithIdx(r).IsInRing():
            continue
        for h in graph[r]:
            atoms = _arm(graph, r, h)
            if h in atoms or any(c in atoms for c in centers if c != r and c not in graph[r]):
                continue
            candidates.setdefault(_key(mol, atoms, r), []).append((r, h, atoms))
    for arms in _rank_candidates(candidates, None):
        if len(arms) != len(centers):
            continue
        name = _attempt(mol, graph, [], arms, "anion")
        if name is not None:
            return name
    return None


def ester_multiplicative_name(mol, ester_oxygens):
    """Identical ester groups joined by a linker that contains heteroatoms: 'oxydi(ethane-2,1-diyl) diacetate'
    (P-65.6.3.2.2)."""
    if len(Chem.GetMolFrags(mol)) != 1 or len(ester_oxygens) < 2:
        return None
    graph = adjacency(mol)
    candidates = {}
    for r in ester_oxygens:
        for h in graph[r]:
            atoms = _arm(graph, r, h)
            if h in atoms or any(c in atoms for c in ester_oxygens if c != r and c not in graph[r]):
                continue
            candidates.setdefault(_key(mol, atoms, r), []).append((r, h, atoms))
    for arms in _rank_candidates(candidates, None):
        if len(arms) != len(ester_oxygens):
            continue
        name = _attempt(mol, graph, [], arms, "ester")
        if name is not None:
            return name
    return None


def _rank_candidates(candidates, anchors, heads=None, mol=None, strict=False):
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
                    if strict:
                        if len(heads & covered) <= _most_groups_on_one_parent(mol, heads):
                            continue
                    elif not anchors <= covered and not _beats_substitutive(mol, heads, covered):
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


def _most_groups_on_one_parent(mol, heads):
    """The most principal groups one acyclic carbon chain or one ring system carries as suffixes; a parent that
    carries as many as the multiplied names express is the parent (P-44.1.1, P-62.2.4.1.3)."""
    graph = adjacency(mol)

    def weight(atom):
        return sum(1 for n in graph[atom] if n in heads)

    best = max(
        (sum(weight(a) for a in system) for system in (_ring_system_of(mol, a.GetIdx()) for a in mol.GetAtoms() if a.IsInRing()) if system),
        default=0,
    )
    chain = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}

    def downward(atom, parent):
        nonlocal best
        legs = sorted((downward(n, atom) for n in graph[atom] if n in chain and n != parent), reverse=True)
        best = max(best, weight(atom) + sum(legs[:2]))
        return weight(atom) + (legs[0] if legs else 0)

    seen = set()
    for start in sorted(chain):
        if start not in seen:
            stack, component = [start], {start}
            while stack:
                for n in graph[stack.pop()]:
                    if n in chain and n not in component:
                        component.add(n)
                        stack.append(n)
            seen |= component
            downward(start, None)
    return best


def _chain_depth(mol, graph, atom, parent):
    best = 0
    for n in graph[atom]:
        neighbor = mol.GetAtomWithIdx(n)
        if n != parent and neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
            best = max(best, _chain_depth(mol, graph, n, atom))
    return 1 + best


def _side_rank(mol, graph, nitrogen, carbon):
    """(ring bound to nitrogen, carbons of the longest chain through `carbon`): the order in which P-44.1.2.2 and
    P-44.3 choose the amine parent among the groups on the nitrogen."""
    if mol.GetAtomWithIdx(carbon).IsInRing():
        return (1, 0)
    depths = sorted(
        (
            _chain_depth(mol, graph, n, carbon)
            for n in graph[carbon]
            if n != nitrogen and mol.GetAtomWithIdx(n).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(n).IsInRing()
        ),
        reverse=True,
    )
    return (0, 1 + sum(depths[:2]))


def _amine_nitrogen_arm(mol, graph, nitrogen, linked):
    """The amine parent held by `nitrogen` is multiplied through that nitrogen when the groups it keeps outrank the
    acyclic carbon `linked` that joins it to the linking group (P-62.2.5.1, P-44.3)."""
    atom = mol.GetAtomWithIdx(nitrogen)
    if (
        atom.IsInRing()
        or atom.GetFormalCharge()
        or atom.GetIsotope()
        or atom.GetNumRadicalElectrons()
        or atom.GetDegree() < 2
        or any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds())
    ):
        return False
    other = mol.GetAtomWithIdx(linked)
    if other.GetAtomicNum() != 6 or other.IsInRing() or is_functional_carbon(mol, linked):
        return False
    kept = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != linked]
    if any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 or is_functional_carbon(mol, n) for n in kept):
        return False
    return max(_side_rank(mol, graph, nitrogen, n) for n in kept) > _side_rank(mol, graph, nitrogen, linked)


def _bare_chain(graph, root, senior_atoms, arm_atoms):
    if arm_atoms != senior_atoms:
        return False
    end = [a for a in senior_atoms if sum(n in senior_atoms for n in graph[a]) <= 1]
    return root in end and all(sum(n in senior_atoms for n in graph[a]) <= 2 for a in senior_atoms)


def _senior_chain(mol, graph, root, senior_atoms, arm_atoms):
    """The linear run of the senior element from `root` when it holds every such atom of the arm, else None."""
    chain = [root]
    while True:
        onward = [n for n in graph[chain[-1]] if n in arm_atoms and n in senior_atoms and n not in chain]
        if len(onward) > 1:
            return None
        if not onward:
            return chain if set(chain) == set(senior_atoms) else None
        chain.append(onward[0])


_HYDRIDE_STEMS = {14: "silane", 32: "germane", 50: "stannane", 82: "plumbane", 15: "phosphane", 33: "arsane"}
_HYDRIDE_VALENCE = {14: 4, 32: 4, 50: 4, 82: 4, 15: 3, 33: 3}


def _substituted_chain_unit(mol, graph, root, arm_atoms):
    """P-15.3.2.4.1: the unit of a multiplicative name that is a substituted dinuclear or longer hydride attached at
    its first atom, with the prefixes of the unit cited inside the multiplied group ('pentafluorodisilane',
    '2-methyldisilane'); None unless it is a linear chain of one element carrying single-bonded groups."""
    from ._common import group_substituents
    from ._numerals import multiplying_prefix
    from ._substituents import format_substituent_prefixes

    atom = mol.GetAtomWithIdx(root)
    z = atom.GetAtomicNum()
    senior = {a for a in arm_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == z}
    chain = _senior_chain(mol, graph, root, senior, arm_atoms)
    if chain is None or len(chain) < 2:
        return None
    if any(
        mol.GetAtomWithIdx(a).GetFormalCharge()
        or mol.GetAtomWithIdx(a).GetIsotope()
        or mol.GetAtomWithIdx(a).GetTotalValence() != _HYDRIDE_VALENCE[z]
        for a in chain
    ):
        return None
    halogens = halogen_substituents(mol)
    entries = {}
    for position, member in enumerate(chain, start=1):
        for n in graph[member]:
            if n in chain or n not in arm_atoms:
                continue
            if mol.GetBondBetweenAtoms(member, n).GetBondTypeAsDouble() != 1.0:
                return None
            entries.setdefault(position, []).append(
                name_branch(graph, n, member, halogens, frozenset(), mol=mol, unsaturated=True)
            )
    if not entries:
        return None
    capacity = _HYDRIDE_VALENCE[z] * len(chain) - 2 * (len(chain) - 1) - 1
    count = sum(len(v) for v in entries.values())
    prefix = format_substituent_prefixes(group_substituents(entries), omit_locants=count == capacity)
    return prefix + multiplying_prefix(len(chain)) + _HYDRIDE_STEMS[z]


def _hydride_candidates(mol, graph):
    if any(a.GetAtomicNum() == 7 for a in mol.GetAtoms()):
        return {}
    present = {a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() in _HYDRIDE_ORDER}
    if any(
        n.GetAtomicNum() in (8, 16, 34, 52) and n.GetDegree() == 1 and n.GetTotalNumHs() == 1
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in _HYDRIDE_ORDER
        for n in a.GetNeighbors()
    ):
        return {}
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
            if h in atoms or (
                inside != {r}
                and not _bare_chain(graph, r, inside, atoms)
                and _senior_chain(mol, graph, r, inside, atoms) is None
            ):
                continue
            candidates.setdefault(_key(mol, atoms, r), []).append((r, h, atoms))
    return candidates


def _amine_candidates(mol, graph):
    """Secondary or tertiary amines whose nitrogen is joined to a linking heteroatom (hydroxylamine O-substitution,
    P-68.3.1.1.1.3): the amine parents are multiplied through nitrogen."""
    candidates = {}
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() != 7
            or atom.IsInRing()
            or atom.GetFormalCharge()
            or atom.GetIsotope()
            or atom.GetNumRadicalElectrons()
            or any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds())
        ):
            continue
        neighbors = list(atom.GetNeighbors())
        joined = [n for n in neighbors if n.GetAtomicNum() != 6]
        if len(joined) != 1 or len(joined) == len(neighbors):
            continue
        linker = joined[0]
        if linker.GetAtomicNum() == 7 or not _is_linker_atom(mol, linker):
            continue
        atoms = _arm(graph, atom.GetIdx(), linker.GetIdx())
        if linker.GetIdx() not in atoms:
            candidates.setdefault(_key(mol, atoms, atom.GetIdx()), []).append((atom.GetIdx(), linker.GetIdx(), atoms))
    return candidates


def _polyamine_candidates(mol, graph):
    """Diamine parents joined through a nitrogen atom by a divalent linking group (P-62.2.4.1.3, P-62.2.5.2): each
    arm is the part of the molecule on the nitrogen side of one N-C bond."""
    candidates = {}
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() != 7
            or atom.IsInRing()
            or atom.GetFormalCharge()
            or atom.GetIsotope()
            or atom.GetDegree() != 2
            or atom.GetTotalNumHs() != 1
            or any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds())
        ):
            continue
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6 or neighbor.IsInRing() or is_functional_carbon(mol, neighbor.GetIdx()):
                continue
            atoms = _arm(graph, atom.GetIdx(), neighbor.GetIdx())
            if neighbor.GetIdx() in atoms or sum(mol.GetAtomWithIdx(a).GetAtomicNum() == 7 for a in atoms) < 2:
                continue
            candidates.setdefault(_key(mol, atoms, atom.GetIdx()), []).append((atom.GetIdx(), neighbor.GetIdx(), atoms))
    return candidates


def _phosphonous_acid_atoms(mol):
    """P, As or Sb atoms of a phosphonous-type group: one carbon neighbour and two hydroxyl groups."""
    return [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in (15, 33, 51)
        and a.GetDegree() == 3
        and not a.GetFormalCharge()
        and sum(
            1
            for n in a.GetNeighbors()
            if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1
            and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        )
        == 2
    ]


def _phosphonic_candidates(mol, graph):
    """Phosphonic, arsonic and stibonic acid groups joined to carbon: the acid outranks amines and ethers (P-41, P-45.1.2)."""
    from ._phosphonic_acid import _SENIOR_ACIDS, _phosphonic_acid_phosphorus_atoms

    if any(mol.HasSubstructMatch(query) for query in _SENIOR_ACIDS):
        return {}
    candidates = {}
    for center in [*_phosphonic_acid_phosphorus_atoms(mol), *_phosphonous_acid_atoms(mol)]:
        oxygens = {n.GetIdx() for n in center.GetNeighbors() if n.GetAtomicNum() == 8}
        roots = [n.GetIdx() for n in center.GetNeighbors() if n.GetIdx() not in oxygens]
        if len(roots) != 1:
            continue
        root = roots[0]
        if mol.GetAtomWithIdx(root).GetAtomicNum() != 6:
            return {}
        atoms = oxygens | {center.GetIdx()}
        candidates.setdefault(_key(mol, atoms, center.GetIdx()), []).append((center.GetIdx(), root, atoms))
    return candidates


_CENTRE_ACID = re.compile(r"(?:bor|phosph|ars|stib)(?:on|in)(?:ic|ous) acid$")
_INORGANIC_AMIDE = re.compile(r"(?:nitramide|(?:ous|ic) amide)$")
_CARBAMIC_ACID = re.compile(r"carbam\w*ic (?:[A-Za-z,]+-)?acid$|carbam\w*ic acid$")
_UREA_STEMS = {8: "urea", 16: "thiourea", 34: "selenourea", 52: "tellurourea"}


def _terminal(atom, symbol, hydrogens=None):
    return (
        atom.GetAtomicNum() == symbol
        and atom.GetDegree() == 1
        and not atom.GetFormalCharge()
        and (hydrogens is None or atom.GetTotalNumHs() == hydrogens)
    )


def _unit_arms(mol, graph, attach, roots):
    """(attach, root, unit atoms) for each carbon `root` that could be the linking atom of the unit at `attach`."""
    found = []
    for root in roots:
        atoms = _arm(graph, attach, root)
        if root not in atoms:
            found.append((attach, root, atoms))
    return found


def _acid_centre_candidates(mol, graph):
    """Boronic, borinic, phosphinic and phosphinous acid groups (and the arsenic and antimony analogues) joined through
    carbon: the centre keeps its hydroxy group(s) and one carbon substituent in the unit (P-67.1.1, P-68.1.4)."""
    from ._phosphonic_acid import _SENIOR_ACIDS

    if any(mol.HasSubstructMatch(query) for query in _SENIOR_ACIDS):
        return {}
    candidates = {}
    for center in mol.GetAtoms():
        if center.GetAtomicNum() not in (5, 15, 33, 51) or center.GetFormalCharge() or center.IsInRing():
            continue
        neighbors = list(center.GetNeighbors())

        def bonded(order, group):
            return [n for n in neighbors if group(n) and mol.GetBondBetweenAtoms(center.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == order]

        hydroxyl = bonded(1.0, lambda n: _terminal(n, 8, 1))
        oxo = bonded(2.0, lambda n: _terminal(n, 8, 0))
        carbons = bonded(1.0, lambda n: n.GetAtomicNum() == 6)
        if len(hydroxyl) + len(oxo) + len(carbons) != len(neighbors):
            continue
        if center.GetAtomicNum() == 5:
            fits = not oxo and hydroxyl and len(hydroxyl) + len(carbons) == 3
        else:
            fits = len(carbons) == 2 and len(hydroxyl) == 1 and len(oxo) <= 1
        if not fits:
            continue
        for found in _unit_arms(mol, graph, center.GetIdx(), [n.GetIdx() for n in carbons]):
            candidates.setdefault(_key(mol, found[2], found[0]), []).append(found)
    return candidates


def _carbamic_core(mol, nitrogen):
    """The carbon of a carbamic or carbamodithioic acid group on `nitrogen`, with its chalcogen atoms, or None."""
    for carbon in nitrogen.GetNeighbors():
        if carbon.GetAtomicNum() != 6 or carbon.GetFormalCharge() or carbon.GetDegree() != 3:
            continue
        doubly = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() in (8, 16) and _terminal(n, n.GetAtomicNum(), 0)
                  and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
        singly = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() in (8, 16) and _terminal(n, n.GetAtomicNum(), 1)
                  and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(doubly) == 1 and len(singly) == 1:
            return {carbon.GetIdx(), doubly[0].GetIdx(), singly[0].GetIdx()}
    return None


def _urea_core(mol, nitrogen):
    """{carbon, chalcogen, far nitrogen} of a urea or thiourea group whose near nitrogen is `nitrogen`, or None."""
    for carbon in nitrogen.GetNeighbors():
        if carbon.GetAtomicNum() != 6 or carbon.GetFormalCharge() or carbon.GetDegree() != 3 or carbon.IsInRing():
            continue
        chalcogen = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() in _UREA_STEMS and _terminal(n, n.GetAtomicNum(), 0)
                     and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
        far = [n for n in carbon.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx() and n.GetAtomicNum() == 7 and not n.GetFormalCharge()
               and not n.IsInRing() and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(chalcogen) == 1 and len(far) == 1:
            return {carbon.GetIdx(), chalcogen[0].GetIdx(), far[0].GetIdx()}
    return None


def _nitrogen_unit_candidates(mol, graph, core_of, senior):
    """Groups of the form R-N(-R')-core joined through the carbon `root` of the nitrogen: `core_of` finds the core atoms
    of the group, and `senior` tells whether a class that outranks the group is present."""
    candidates, cores = {}, set()
    for nitrogen in mol.GetAtoms():
        if nitrogen.GetAtomicNum() != 7 or nitrogen.IsInRing() or nitrogen.GetFormalCharge() or nitrogen.GetIsotope():
            continue
        core = core_of(mol, nitrogen)
        if core is None:
            continue
        cores |= core | {nitrogen.GetIdx()}
        roots = [
            n.GetIdx()
            for n in nitrogen.GetNeighbors()
            if n.GetAtomicNum() == 6 and n.GetIdx() not in core and not is_functional_carbon(mol, n.GetIdx())
            and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        for found in _unit_arms(mol, graph, nitrogen.GetIdx(), roots):
            candidates.setdefault(_key(mol, found[2], found[0]), []).append(found)
    return {} if senior(mol, cores) else candidates


def _inorganic_amide_core(mol, nitrogen):
    """The atoms of the nitroso, nitro or halogen group that makes `nitrogen` the amide nitrogen of nitrous, nitric or a
    halogen oxoacid, or None when it carries no such group or several."""
    from ._halogen_amide import _oxo_oxygens
    from ._nitro_group import nitro_or_nitroso_group

    cores = []
    for neighbor in nitrogen.GetNeighbors():
        if neighbor.GetAtomicNum() in HALOGEN_PREFIXES and not neighbor.GetFormalCharge():
            oxygens = _oxo_oxygens(neighbor, nitrogen)
            if oxygens is not None:
                cores.append({neighbor.GetIdx(), *(o.GetIdx() for o in oxygens)})
            continue
        group = nitro_or_nitroso_group(mol, neighbor)
        if group is not None and group[0] in ("nitro", "nitroso"):
            cores.append({neighbor.GetIdx(), *(o.GetIdx() for o in group[1])})
    return cores[0] if len(cores) == 1 else None


def _senior_to_inorganic_amide(mol, cores):
    from ._inorganic_acid_derivative import _only_prefix_groups_beyond

    return not _only_prefix_groups_beyond(mol, cores)


def _senior_to_carbamic(mol, cores):
    from ._chalcogenourea import _SENIOR_TO_UREA

    return any(match[0] not in cores for query in (_SENIOR_TO_UREA[0], _SENIOR_TO_UREA[-1]) for match in mol.GetSubstructMatches(query))


def _senior_to_urea(mol, cores):
    from ._chalcogenourea import has_group_senior_to_urea

    return has_group_senior_to_urea(mol, cores)


def _urea_unit_text(mol, graph, attach, root, atoms):
    """'urea', 'N'-methylurea', 'thiourea': the unit joined at the near nitrogen `attach`, which takes the locant N
    (P-66.1.6.1.1.2); the far nitrogen takes N'."""
    from ._common import group_substituents
    from ._substituents import format_substituent_prefixes

    core = _urea_core(mol, mol.GetAtomWithIdx(attach))
    if core is None:
        return None
    carbon = next(a for a in core if mol.GetAtomWithIdx(a).GetAtomicNum() == 6)
    chalcogen = next(a for a in core if mol.GetAtomWithIdx(a).GetAtomicNum() in _UREA_STEMS)
    far = next(a for a in core if a not in (carbon, chalcogen))
    halogens = halogen_substituents(mol)
    entries = {}
    for locant, nitrogen, excluded in (("N", attach, (root, carbon)), ("N'", far, (carbon,))):
        for n in graph[nitrogen]:
            if n not in excluded:
                entries.setdefault(locant, []).append(name_branch(graph, n, nitrogen, halogens, mol=mol, unsaturated=True))
    prefix = format_substituent_prefixes(group_substituents(entries)) if entries else ""
    return prefix + _UREA_STEMS[mol.GetAtomWithIdx(chalcogen).GetAtomicNum()]


def chain_multiplicative_name(mol, stereo):
    name = _chain_multiplicative_name(mol, stereo)
    return name.replace("(1,2-dioxoethane-1,2-diyl)", "oxalyl") if name else name


def _chain_multiplicative_name(mol, stereo):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    graph = adjacency(mol)
    for arms in _rank_candidates(_phosphonic_candidates(mol, graph), None):
        name = _attempt(mol, graph, stereo, arms, "phosphonic")
        if name is not None:
            return name
    for kind, candidates in (
        ("centre_acid", _acid_centre_candidates(mol, graph)),
        ("carbamic", _nitrogen_unit_candidates(mol, graph, _carbamic_core, _senior_to_carbamic)),
        ("urea", _nitrogen_unit_candidates(mol, graph, _urea_core, _senior_to_urea)),
        ("inorganic_amide", _nitrogen_unit_candidates(mol, graph, _inorganic_amide_core, _senior_to_inorganic_amide)),
    ):
        for arms in _rank_candidates(candidates, None):
            name = _attempt(mol, graph, stereo, arms, kind)
            if name is not None:
                return name
    found = _principal_atoms(mol)
    heterone_oxygens = chalcogen_chain_oxygen_count(mol)
    if heterone_oxygens:
        # P-64.7.2: a heterone ranks with the ketones and the one holding more doubly bonded oxygen atoms is the parent
        from ._polyfunctional import _SENIORITY

        if found is None or _SENIORITY.index(found[0]) > _SENIORITY.index("ketone"):
            return None
    if found is None:
        for arms in _rank_candidates(_amine_candidates(mol, graph), None):
            name = _attempt(mol, graph, stereo, arms, "amine")
            if name is not None:
                return name
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
    from ._polyfunctional import _SENIORITY

    acyl_prefixes = (
        frozenset(
            a.GetIdx()
            for a in mol.GetAtoms()
            if _acyl_carbon(a)
            and not any(
                _acyl_branch(graph, a.GetIdx(), n.GetIdx()) & anchors for n in a.GetNeighbors() if n.GetAtomicNum() == 7
            )
        )
        if _SENIORITY.index(principal) < _SENIORITY.index("amide")
        else frozenset()
    )
    linkers = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if _is_linker_atom(mol, a, acyl_prefixes) and not (principal == "amine" and a.GetAtomicNum() == 7)
    }
    candidates, nitrogen, amine_nitrogen = {}, {}, {}
    for r in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(r)
        if atom.IsInRing():
            continue
        for h in graph[r]:
            bond = mol.GetBondBetweenAtoms(r, h)
            if atom.GetAtomicNum() == 6:
                ok = (h in linkers and not (principal == "ketone" and r in anchors)) or (
                    not (principal == "aldehyde" and r in anchors)
                    and mol.GetAtomWithIdx(h).IsInRing()
                    and bond.GetBondTypeAsDouble() in (1.0, 2.0)
                )
                target = candidates
            elif atom.GetAtomicNum() == 7 and principal in ("amide", "sulfonamide", "imine") and r in anchors:
                other = mol.GetAtomWithIdx(h)
                ok = bond.GetBondTypeAsDouble() == 1.0 and (
                    (other.GetAtomicNum() == 6 and not is_functional_carbon(mol, h))
                    or (other.GetAtomicNum() not in (6, 7) and h in linkers)
                )
                target = nitrogen
            elif atom.GetAtomicNum() == 7 and principal == "amine" and r in anchors:
                ok = _amine_nitrogen_arm(mol, graph, r, h)
                target = amine_nitrogen
            else:
                continue
            if ok:
                atoms = _arm(graph, r, h)
                if h not in atoms:
                    target.setdefault(_key(mol, atoms, r), []).append((r, h, atoms))
    for arms in _rank_candidates(candidates, anchors, heads, mol):
        if principal == "ketone" and heterone_oxygens > max(len(heads & atoms) for _, _, atoms in arms):
            continue
        name = _attempt(mol, graph, stereo, arms)
        if name is not None:
            return name
    for arms in _rank_candidates(nitrogen, anchors, heads, mol):
        name = _attempt(mol, graph, stereo, arms, "imine" if principal == "imine" else "amide")
        if name is not None:
            return name
    for arms in _rank_candidates(amine_nitrogen, anchors, heads, mol, strict=True):
        name = _attempt(mol, graph, stereo, arms, "amine")
        if name is not None:
            return name
    if principal == "amine":
        for arms in _rank_candidates(_polyamine_candidates(mol, graph), anchors, heads, mol):
            name = _attempt(mol, graph, stereo, arms, "polyamine")
            if name is not None:
                return name
    return None


def has_chain_multiplicative_shape(mol):
    try:
        return chain_multiplicative_name(mol, []) is not None
    except Exception:
        return False
