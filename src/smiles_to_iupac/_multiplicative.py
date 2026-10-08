"""Multiplicative names for identical ring parent structures joined by a di- or
polyvalent group (P-15.3, P-51.3): '1,1'-methylenedibenzene',
'4,4'-oxydibenzoic acid', '1,1',1''-(ethane-1,1,2-triyl)tribenzene'.
The units must be the senior parent (P-41 class, then more numerous), attached
by identical single bonds to one central group, optionally via identical arms.
"""

import itertools
import re
from dataclasses import dataclass

from rdkit import Chem

from ._common import UnsupportedStructure, specified_stereo_elements
from ._multiplicative_groups import SUFFIX_RANKS, classify
from ._ring_system_seniority import ring_seniority_key
from ._multiplicative_linker import DecompositionRejected, name_component
from ._multiplicative_prefix import SUFFIX_CARRIERS, hook_suspended
from ._multiplicative_ring import (
    bare_polycyclic_unit,
    name_monocyclic_unit,
    principal_rank_of,
    spec_of,
    substituted_polycyclic_unit,
)
from ._multiplicative_text import enclose, multiplier_word, primed_locants

_MAX_GROUP = 8
_SENIOR_HYDRIDE_ATOMS = {5, 14, 15, 32, 33, 50, 51, 82, 83}


@dataclass
class _Context:
    groups: list
    suffix_group: object
    name_function: object
    stereo: list
    used: set
    entry_cache: dict
    entry: object = None
    soft: bool = False


@dataclass
class _Unit:
    node: tuple
    ring_atoms: frozenset
    atoms: frozenset
    junction: int
    linker_atom: int
    key: str


def _ring_systems(mol):
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ring in mol.GetRingInfo().AtomRings():
        for a in ring:
            parent.setdefault(a, a)
        for a in ring[1:]:
            parent[find(a)] = find(ring[0])
    systems = {}
    for a in parent:
        systems.setdefault(find(a), set()).add(a)
    return list(systems.values())


def _fragment_key(mol, atoms, marked):
    rw = Chem.RWMol(mol)
    atom = rw.GetAtomWithIdx(marked)
    atom.SetAtomMapNum(1)
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(0)
    rw.UpdatePropertyCache(strict=False)
    return Chem.MolFragmentToSmiles(rw, atomsToUse=sorted(atoms), canonical=True)


def _bare_key(mol, atoms):
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
    return Chem.MolFragmentToSmiles(rw, atomsToUse=sorted(atoms), canonical=True)


def _build_tree(mol):
    systems = _ring_systems(mol)
    node_of = {}
    for i, atoms in enumerate(systems):
        for a in atoms:
            node_of[a] = ("R", i)
    for atom in mol.GetAtoms():
        node_of.setdefault(atom.GetIdx(), ("A", atom.GetIdx()))
    adj = {}
    for bond in mol.GetBonds():
        x, y = node_of[bond.GetBeginAtomIdx()], node_of[bond.GetEndAtomIdx()]
        if x != y:
            adj.setdefault(x, set()).add(y)
            adj.setdefault(y, set()).add(x)
    for n in set(node_of.values()):
        adj.setdefault(n, set())
    edges = sum(len(v) for v in adj.values()) // 2
    if edges != len(adj) - 1:
        return None
    return systems, node_of, adj


def _prune(adj, keep):
    remaining = {n: set(v) for n, v in adj.items()}
    queue = [n for n, v in remaining.items() if len(v) <= 1 and n not in keep]
    while queue:
        n = queue.pop()
        if n not in remaining:
            continue
        for m in remaining[n]:
            remaining[m].discard(n)
            if len(remaining[m]) <= 1 and m not in keep:
                queue.append(m)
        del remaining[n]
    return remaining


def _peel_acyclic(adj):
    remaining = {n: set(v) for n, v in adj.items()}
    queue = [n for n, v in remaining.items() if n[0] == "A" and len(v) <= 1]
    while queue:
        n = queue.pop()
        if n not in remaining:
            continue
        for m in remaining[n]:
            remaining[m].discard(n)
            if m[0] == "A" and len(remaining[m]) <= 1:
                queue.append(m)
        del remaining[n]
    return remaining


def _leaf_units(mol, systems, node_of, core):
    units = []
    for node, neighbors in core.items():
        if node[0] != "R" or len(neighbors) != 1:
            continue
        (other,) = neighbors
        ring_atoms = frozenset(systems[node[1]])
        if other[0] == "A":
            linker_atom = other[1]
            junction = next(
                (n.GetIdx() for n in mol.GetAtomWithIdx(linker_atom).GetNeighbors() if n.GetIdx() in ring_atoms),
                None,
            )
        else:
            neighbor_atoms = systems[other[1]]
            bridge = next(
                (
                    (a, n.GetIdx())
                    for a in ring_atoms
                    for n in mol.GetAtomWithIdx(a).GetNeighbors()
                    if n.GetIdx() in neighbor_atoms
                ),
                None,
            )
            if bridge is None or _bare_key(mol, ring_atoms) == _bare_key(mol, neighbor_atoms):
                continue
            junction, linker_atom = bridge
        if junction is None:
            continue
        bond = mol.GetBondBetweenAtoms(junction, linker_atom)
        order = bond.GetBondTypeAsDouble()
        if order == 2 and mol.GetAtomWithIdx(junction).GetIsAromatic():
            continue
        if order not in (1, 2):
            continue
        atoms = set(ring_atoms)
        stack = list(ring_atoms)
        while stack:
            a = stack.pop()
            for n in mol.GetAtomWithIdx(a).GetNeighbors():
                idx = n.GetIdx()
                if idx in atoms or node_of[idx] in core:
                    continue
                atoms.add(idx)
                stack.append(idx)
        units.append(_Unit(node, ring_atoms, frozenset(atoms), junction, linker_atom, _fragment_key(mol, atoms, junction)))
    return units


def _class_gate(mol, groups, selected, systems, node_of):
    unit_atoms = set().union(*(u.atoms for u in selected))
    best_in = min((g.rank for g in groups if g.anchor in unit_atoms and g.rank <= SUFFIX_RANKS["amine"]), default=None)
    best_out = min((g.rank for g in groups if g.anchor not in unit_atoms and g.rank <= SUFFIX_RANKS["amine"]), default=None)
    if best_in is not None:
        if best_out is None or best_out > best_in:
            return True
        # P-44.1.1: at equal class the parent holding more principal groups wins; a tie goes to the rings (P-44.1.2.2)
        inside = sum(1 for g in groups if g.anchor in unit_atoms and g.rank == best_in)
        outside = sum(1 for g in groups if g.anchor not in unit_atoms and g.rank == best_in)
        return best_out == best_in and inside >= outside
    if best_out is not None:
        return False
    selected_nodes = {u.node for u in selected}
    unit_key = ring_seniority_key(mol, selected[0].ring_atoms)
    return all(
        ring_seniority_key(mol, atoms) >= unit_key for i, atoms in enumerate(systems) if ("R", i) not in selected_nodes
    )


_HOMONUCLEAR_RUNS = (7, 8, 14, 15, 16, 32, 33, 34, 50, 51, 52, 82, 83)


def _components(mol, linker_nodes, systems):
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    atoms = [n[1] for n in linker_nodes if n[0] == "A"]
    for a in atoms:
        parent[a] = a
    atom_set = set(atoms)
    for a in atoms:
        za = mol.GetAtomWithIdx(a).GetAtomicNum()
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            b = n.GetIdx()
            if b not in atom_set:
                continue
            zb = n.GetAtomicNum()
            order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
            single = order == 1
            if (za == 6 and zb == 6) or (za == zb and za in _HOMONUCLEAR_RUNS and (single or (za == 7 and order == 2))):
                parent[find(a)] = find(b)
    comps = {}
    for a in atoms:
        comps.setdefault(("C", find(a)), []).append(a)
    result = {}
    for cid, members in comps.items():
        z = mol.GetAtomWithIdx(members[0]).GetAtomicNum()
        result[cid] = ("carbon" if z == 6 else "hetero", sorted(members))
    for n in linker_nodes:
        if n[0] == "R":
            result[("R", n[1])] = ("ring", sorted(systems[n[1]]))
    return result


def _branch_atoms(mol, cut_a, cut_b):
    seen = {cut_b}
    stack = [cut_b]
    while stack:
        a = stack.pop()
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            idx = n.GetIdx()
            if idx in seen or (a == cut_b and idx == cut_a):
                continue
            seen.add(idx)
            stack.append(idx)
    return seen


def _comp_of_atom(components):
    result = {}
    for cid, (_, atoms) in components.items():
        for a in atoms:
            result[a] = cid
    return result


def _find_center(mol, components, selected):
    comp_of = _comp_of_atom(components)
    junctions = {u.junction: i for i, u in enumerate(selected)}
    edges = {cid: [] for cid in components}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = int(bond.GetBondTypeAsDouble())
        for x, y in ((a, b), (b, a)):
            if x in comp_of and y in comp_of and comp_of[x] != comp_of[y]:
                edges[comp_of[x]].append((x, y, order))
            elif x in comp_of and y in junctions:
                edges[comp_of[x]].append((x, y, order))
    valid = []
    for cid, cedges in edges.items():
        keys = set()
        ok = True
        for x, y, order in cedges:
            branch = _branch_atoms(mol, x, y)
            if sum(1 for j in junctions if j in branch) != 1:
                ok = False
                break
            keys.add((0 if y in junctions else order, _fragment_key(mol, branch, y)))
        if ok and len(keys) == 1 and len(cedges) >= 2:
            valid.append(cid)
    if not valid:
        return None
    return valid, edges, comp_of


def _arm_chain(mol, center, edges, comp_of, first_edge, junction_atoms):
    x, y, order = first_edge
    chain = []
    prev_atom_in_prev, current_atom = x, y
    while current_atom not in junction_atoms:
        cid = comp_of[current_atom]
        cedges = edges[cid]
        toward_center = (current_atom, prev_atom_in_prev, order)
        onward = [e for e in cedges if not (e[0] == current_atom and e[1] == prev_atom_in_prev)]
        if len(onward) != 1:
            return None
        ox, oy, oorder = onward[0]
        chain.append((cid, toward_center, (ox, oy, oorder)))
        prev_atom_in_prev, current_atom, order = ox, oy, oorder
    return chain


def _unit_text(mol, unit, groups, name_function):
    ring_count = sum(1 for ring in mol.GetRingInfo().AtomRings() if set(ring) <= unit.ring_atoms)
    if ring_count == 1 and unit.atoms == unit.ring_atoms and spec_of(mol, unit.ring_atoms) is None:
        return bare_polycyclic_unit(mol, unit.ring_atoms, unit.junction, name_function)
    if ring_count == 1:
        return name_monocyclic_unit(mol, unit.ring_atoms, unit.junction, unit.linker_atom, groups, unit.atoms, name_function)
    if unit.atoms != unit.ring_atoms:
        return substituted_polycyclic_unit(mol, unit.ring_atoms, unit.atoms, unit.junction, name_function)
    return bare_polycyclic_unit(mol, unit.ring_atoms, unit.junction, name_function)


def _attempt(mol, groups, selected, tree, core, name_function):
    systems, node_of, adj = tree
    selected_nodes = {u.node for u in selected}
    span = _prune(core, selected_nodes)
    linker_nodes = [n for n in span if n not in selected_nodes]
    if not linker_nodes or not _class_gate(mol, groups, selected, systems, node_of):
        return None
    components = _components(mol, linker_nodes, systems)
    found = _find_center(mol, components, selected)
    if found is None:
        return None
    valid, edges, comp_of = found
    center = valid[0]

    unit_atoms = set().union(*(u.atoms for u in selected))
    principal = principal_rank_of(groups, unit_atoms)
    if principal is None and any(a.GetAtomicNum() in _SENIOR_HYDRIDE_ATOMS for a in mol.GetAtoms()):
        return None
    if principal is None and any(g.name == "imine" for g in groups):
        return None
    if principal is None and any(
        len(members) >= 2 and mol.GetAtomWithIdx(members[0]).GetAtomicNum() == 7
        for component_kind, members in components.values()
        if component_kind == "hetero"
    ):
        return None
    principal_group = None
    if principal is not None:
        principal_group = next(g.name for g in groups if g.anchor in unit_atoms and g.rank == principal)
        if principal_group not in SUFFIX_CARRIERS and principal_group != "ketone":
            return None
    if principal_group == "amine" and any(mol.GetAtomWithIdx(u.linker_atom).GetAtomicNum() == 7 for u in selected):
        return None
    soft = any(mol.GetBondBetweenAtoms(u.junction, u.linker_atom).GetBondTypeAsDouble() != 1 for u in selected)
    ctx = _Context(groups, principal_group, name_function, specified_stereo_elements(mol) or [], set(), {}, soft=soft)
    stereo = ctx.stereo

    unit = _unit_text(mol, selected[0], groups, name_function)
    if unit is None:
        return None

    kind, atoms = components[center]
    center_edges = edges[center]
    junction_atoms = {u.junction for u in selected}
    arm_chain = _arm_chain(mol, center, edges, comp_of, center_edges[0], junction_atoms)
    if arm_chain is None:
        return None
    try:
        central = name_component(mol, kind, atoms, center_edges, ctx)
        arm_parts = []
        for cid, toward_center, toward_unit in arm_chain:
            akind, aatoms = components[cid]
            arm_parts.append(
                name_component(
                    mol, akind, aatoms, [toward_center, toward_unit], ctx, directed=(toward_unit[0], toward_center[0])
                )
            )
    except DecompositionRejected:
        return None
    except UnsupportedStructure:
        if soft:
            return None
        raise

    represented = set(atoms)
    for cid, _, _ in arm_chain:
        represented.update(components[cid][1])
    elsewhere = set().union(*(set(a) for _, a in components.values())) - represented
    for element_kind, idx, _ in stereo:
        if element_kind == "atom":
            involved = {idx}
        else:
            bond = mol.GetBondWithIdx(idx)
            involved = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
        if (element_kind, idx) not in ctx.used and not involved <= elsewhere:
            raise UnsupportedStructure("stereodescriptors outside the linking group of a multiplicative name are not supported yet")
    return _assemble(len(selected), unit, central, arm_parts)


_MULTIPLIED_HYDRIDE_ARM = re.compile(r"(?:di|tri|tetra|penta|hexa|hepta|octa)(?:silane|germane|stannane|plumbane|phosphane|arsane)")


def _assemble(count, unit, central, arm_parts):
    if not arm_parts:
        linker = enclose(central.text) if (central.has_prefix or central.has_locants) else central.text
    else:
        central_text = enclose(central.text) if central.has_prefix else central.text
        pieces = [enclose(p.text) if p.has_prefix and len(arm_parts) > 1 else p.text for p in arm_parts]
        arm_text = "".join(pieces)
        single_locant_only = (
            len(arm_parts) == 1
            and arm_parts[0].has_locants
            and not arm_parts[0].has_prefix
            and not _MULTIPLIED_HYDRIDE_ARM.match(arm_text)
        )
        multiplier = multiplier_word(count, use_bis=not single_locant_only)
        linker = enclose(central_text + multiplier + enclose(arm_text))
    if unit.substituted or unit.text[0].isdigit():
        unit_text = multiplier_word(count, True) + enclose(unit.text)
    elif unit.has_locants:
        unit_text = multiplier_word(count, False) + enclose(unit.text)
    else:
        unit_text = multiplier_word(count, False) + unit.text
    return f"{primed_locants(unit.junction_locant, count)}-{linker}{unit_text}"


def _is_linear_phane(core):
    rings = sum(1 for n in core if n[0] == "R")
    return rings >= 4 and len(core) >= 7 and all(len(v) <= 2 for v in core.values())


def name_if_multiplicative(mol, name_function=None):
    name = _name_if_multiplicative(mol, name_function)
    return name.replace("(1,2-dioxoethane-1,2-diyl)", "oxalyl") if name else name


def _name_if_multiplicative(mol, name_function=None):
    """The multiplicative name of `mol` when it is one, else None. Raises
    `UnsupportedStructure` when the structure is multiplicative but a
    component of the name can't be built yet."""
    if hook_suspended() or len(Chem.GetMolFrags(mol)) != 1:
        return None
    if mol.GetRingInfo().NumRings() < 2 or any(a.GetFormalCharge() and a.IsInRing() for a in mol.GetAtoms()):
        return None
    tree = _build_tree(mol)
    if tree is None:
        return None
    groups = classify(mol)
    if groups is None:
        return None
    systems, node_of, adj = tree
    core = _peel_acyclic(adj)
    units = _leaf_units(mol, systems, node_of, core)
    by_key = {}
    for u in units:
        by_key.setdefault(u.key, []).append(u)
    candidates = []
    for key, members in by_key.items():
        if len(members) < 2 or len(members) > _MAX_GROUP:
            continue
        rank = ring_seniority_key(mol, members[0].ring_atoms)
        principal = principal_rank_of(groups, members[0].atoms)
        for size in range(len(members), 1, -1):
            for subset in itertools.combinations(members, size):
                candidates.append(((-size, principal if principal is not None else 99, rank, key), list(subset)))
    candidates.sort(key=lambda c: c[0])
    for _, selected in candidates:
        name = _attempt(mol, groups, selected, tree, core, name_function)
        if name is not None:
            if _is_linear_phane(core):
                from ._linear_phane import has_linear_phane_shape, name_linear_phane

                if has_linear_phane_shape(mol):
                    return name_linear_phane(mol)
                raise UnsupportedStructure("the PIN of a linear phane (P-52.2.5.1) is a phane name, not supported yet")
            return name
    return None
