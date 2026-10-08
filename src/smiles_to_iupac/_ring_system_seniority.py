"""Seniority key of a ring system as a parent structure (P-44.2): the general criteria of P-44.2.1, the order of the system
types of P-44.2.2.2 and, within a type, the criteria of P-44.2.2.2.1 to P-44.2.2.2.7. A smaller key is the senior parent.
"""

import re

from rdkit import Chem

from ._common import nonstandard_bonding, UnsupportedStructure

# P-44.2.1 (c) and (g): heteroatoms by seniority, nitrogen taking its place only in (g)
_SENIORITY = (9, 17, 35, 53, 8, 16, 34, 52, 7, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81)
_AHEAD_OF_NITROGEN = {z: i for i, z in enumerate(z for z in _SENIORITY if z != 7)}

# P-44.2.2.2 (a)-(g)
_TYPE_ORDER = ("spiro", "phane", "fused", "bridged_fused", "von_baeyer", "fullerene", "assembly")
_MONOCYCLE = -1


def ring_count(mol, atoms):
    """Rings of the subgraph of ring bonds on `atoms`: bonds - atoms + components."""
    atoms = set(atoms)
    parent = {a: a for a in atoms}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    bonds = 0
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in atoms and b in atoms and bond.IsInRing():
            bonds += 1
            parent[find(a)] = find(b)
    return bonds - len(atoms) + len({find(a) for a in atoms})


def general_key(mol, atoms, as_carbon=frozenset()):
    """P-44.2.1 (a)-(g): heterocycle, nitrogen, earlier heteroatom, rings, skeletal atoms, heteroatoms, earlier heteroatoms;
    the `as_carbon` atoms count as carbon."""
    hetero = [z for z in (mol.GetAtomWithIdx(a).GetAtomicNum() for a in atoms if a not in as_carbon) if z != 6]
    nitrogen = 7 in hetero
    ahead = min((_AHEAD_OF_NITROGEN.get(z, len(_AHEAD_OF_NITROGEN)) for z in hetero if z != 7), default=0)
    return (
        0 if hetero else 1,
        0 if nitrogen else 1,
        0 if nitrogen else ahead,
        -ring_count(mol, atoms),
        -len(atoms),
        -len(hetero),
        tuple(-sum(1 for z in hetero if z == e) for e in _SENIORITY),
    )


def multiple_bond_count(mol, atoms):
    """Double and triple bonds within `atoms`; an aromatic bond counts half where the molecule cannot be kekulized."""
    kekulized = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekulized, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return sum(
            0.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble() > 1
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
        )
    return sum(
        1
        for b in kekulized.GetBonds()
        if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms and b.GetBondTypeAsDouble() > 1
    )


def _hashable(value):
    return tuple(_hashable(v) for v in value) if isinstance(value, (list, tuple)) else value


def ring_seniority_key(mol, ring_atoms, assembly=False):
    """General key, system type and the criteria of that type, then the multiple bonds of P-44.4.1.1; `ring_atoms` of a
    ring assembly hold every component."""
    atoms = frozenset(ring_atoms)
    bonding = sorted((nonstandard_bonding(mol.GetAtomWithIdx(a)) or 0 for a in atoms), reverse=True)
    bonding = [n for n in bonding if n]
    return _parent_key(mol, atoms, assembly) + (-multiple_bond_count(mol, atoms), -len(bonding), tuple(-n for n in bonding))


def _parent_key(mol, atoms, assembly=False, as_carbon=frozenset()):
    """The key of the ring system as a parent hydride, which does not depend on its state of hydrogenation."""
    kind = "assembly" if assembly else system_type(mol, atoms)
    rank = _TYPE_ORDER.index(kind) if kind in _TYPE_ORDER else _MONOCYCLE
    try:
        specific = _TYPE_KEYS[kind](mol, atoms) if kind in _TYPE_KEYS else ()
    except (UnsupportedStructure, ValueError, RuntimeError):
        specific = ()
    return general_key(mol, atoms, as_carbon) + (rank,) + _hashable(specific)


def system_type(mol, atoms):
    """The type of P-44.2.2.2 the ring system `atoms` has as a parent hydride: spiro, phane, fused, bridged_fused,
    von_baeyer or monocycle; a fullerene cage has a type of its own (P-27)."""
    from ._fullerene import is_fullerene_cage
    from ._ring_diyl_numbering import _plain_fused
    from ._spiro_union import _components

    rings = ring_count(mol, atoms)
    if rings == 1:
        return "monocycle"
    if is_fullerene_cage(mol, atoms):
        return "fullerene"
    if _phane_key(mol, atoms) is not None:
        return "phane"
    if sum(1 for c in _components(mol) if c["atoms"] <= atoms) > 1:
        return "spiro"
    # P-52.2.4.1: fusion nomenclature needs two rings of five or more members
    if _plain_fused(mol, atoms) and sum(len(r) >= 5 for r in Chem.GetSymmSSSR(_skeleton(mol, atoms))) >= 2:
        return "fused"
    if rings >= 3 and _bridged_parents(mol, atoms):
        return "bridged_fused"
    return "von_baeyer"


def _phane_key(mol, atoms):
    """P-44.2.2.2.2 key when the ring system holds two disjoint rings and is a parent cyclic phane, else None."""
    from ._phane_general import phane_seniority_key
    from ._ring_diyl_numbering import _bare_skeleton

    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= atoms]
    if not any(not a & b for i, a in enumerate(rings) for b in rings[i + 1 :]):
        return None
    try:
        return phane_seniority_key(_bare_skeleton(mol, atoms)[0])
    except (UnsupportedStructure, ValueError, RuntimeError):
        return None


def _bridged_parents(mol, atoms):
    from ._fusion_bridged import bridged_parents

    try:
        return bridged_parents(mol, atoms)
    except UnsupportedStructure:
        return []


def _skeleton(mol, atoms):
    from ._fusion_components import skeleton

    order = sorted(atoms)
    index = {a: i for i, a in enumerate(order)}
    bonds = [
        (index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()])
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in index and b.GetEndAtomIdx() in index
    ]
    return skeleton([mol.GetAtomWithIdx(a).GetSymbol() for a in order], bonds)


def _component_keys(root):
    keys, stack = [], [root.root if hasattr(root, "root") else root]
    while stack:
        part = stack.pop()
        keys.append(part.comp.senior_key)
        stack.extend(part.children)
    return sorted(keys)


def _fused_key(mol, atoms):
    """P-44.2.2.2.3: larger rings first, more rings in a horizontal row, fusion descriptor letters then numbers, senior
    components in order of decreasing seniority."""
    from ._fused_numbering import orientation_key
    from ._fusion_name import fusion_name_keyed

    sk = _skeleton(mol, atoms)
    sizes = sorted((len(r) for r in Chem.GetSymmSSSR(sk)), reverse=True)
    _, root, descriptor = fusion_name_keyed(sk)
    return (tuple(-s for s in sizes), orientation_key(sk)[0], descriptor, tuple(_component_keys(root)))


def _bridged_key(mol, atoms):
    """P-44.2.2.2.4 (a)-(n) of the best name of the bridged fused system."""
    parents = _bridged_parents(mol, atoms)
    if not parents:
        raise UnsupportedStructure("no bridged fused name")
    return min(_bridged_parent_key(mol, atoms, p) for p in parents)


def _bridged_parent_key(mol, atoms, parent):
    from ._fused_numbering import HETERO_RANK, _locant_key

    fused = frozenset(atoms) - parent.bridge_atoms
    locant = {a: _locant_key(str(loc)) for a, loc in parent.position_of.items()}
    parts = parent.parts
    attachments = [
        sorted(
            locant[n.GetIdx()]
            for a in part.atoms
            for n in mol.GetAtomWithIdx(a).GetNeighbors()
            if n.GetIdx() not in part.atoms and n.GetIdx() in locant
        )
        for part in parts
    ]
    symbols = {a: mol.GetAtomWithIdx(a).GetSymbol() for a in parent.bridge_atoms}
    bridge_hetero = [(locant[a], symbol) for a, symbol in symbols.items() if symbol != "C"]
    by_element = [
        tuple(sorted(loc for loc, s in bridge_hetero if s == element))
        for element in sorted({s for _, s in bridge_hetero}, key=HETERO_RANK.get)
    ]
    independent = [i for i, p in enumerate(parts) if p.independent]
    dependent = [i for i, p in enumerate(parts) if not p.independent]
    return (
        -ring_count(mol, fused),
        -len(fused),
        sum(1 for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() != 6),
        _parent_key(mol, fused),
        tuple(sorted(loc for locs in attachments for loc in locs)),
        tuple(sorted(loc for loc, _ in bridge_hetero)),
        tuple(by_element),
        sum(1 for p in parts if len(p.units) > 1),
        len(dependent),
        sum(len(parts[i].atoms) for i in dependent),
        -sum(1 for p in parts if p.valence == 2),
        tuple(sorted(loc for i in independent for loc in attachments[i])),
        tuple(sorted(loc for i in dependent for loc in attachments[i])),
        -multiple_bond_count(mol, fused),
    )


def _von_baeyer_key(mol, atoms):
    """P-44.2.2.2.5: ring-size descriptor in order of citation, then the superscript locants as a set and in citation order."""
    from ._bicyclic import find_bicyclic_core
    from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates
    from ._ring_diyl_numbering import _bare_skeleton

    bare, _, _ = _bare_skeleton(mol, atoms)
    rings = ring_count(bare, range(bare.GetNumAtoms()))
    if rings == 2:
        core = find_bicyclic_core(bare)
        if core is None:
            raise UnsupportedStructure("no bicyclic core")
        return (tuple(sorted((len(b) for b in core[2]), reverse=True)), (), ())
    core = find_polycyclic_core(bare, rings)
    if core is None:
        raise UnsupportedStructure("no polycyclic core")
    _, parent, _ = min(iter_polycyclic_candidates(core, rings), key=lambda c: c[2])
    descriptor = re.search(r"\[(.*)\]", parent).group(1).split(".")
    sizes = [int(x.split("^")[0]) for x in descriptor]
    superscripts = [int(n) for x in descriptor if "^" in x for n in x.split("^")[1].split(",")]
    return (tuple(sizes), tuple(sorted(superscripts)), tuple(superscripts))


def _replacement_atoms(mol, comp):
    """The heteroatoms of a spiro component that the name introduces by skeletal replacement, so that the component is
    compared as a hydrocarbon (P-24.5.2)."""
    from ._spiro_union import _component

    hetero = frozenset(a for a in comp["atoms"] if mol.GetAtomWithIdx(a).GetAtomicNum() != 6)
    if not hetero:
        return hetero
    try:
        return hetero if _component(mol, comp, False)["replacement"] else frozenset()
    except (UnsupportedStructure, ValueError, RuntimeError):
        return frozenset()


def _spiro_key(mol, atoms):
    """P-44.2.2.2.1: more spiro fusions, saturated monocyclic components with lower spiro locants, then discrete
    components by seniority and in order of citation, with the spiro atoms of lower locants."""
    from ._polyspiro_union import has_spiro_union_shape, spiro_union_numberings
    from ._ring_diyl_numbering import _von_baeyer
    from ._spiro_union import _components

    comps = [c for c in _components(mol) if c["atoms"] <= atoms]
    owners = {a: [i for i, c in enumerate(comps) if a in c["atoms"]] for a in atoms}
    spiro_atoms = {a for a, o in owners.items() if len(o) > 1}
    saturated = all(c["rings"] == 1 for c in comps) and multiple_bond_count(mol, atoms) == 0
    graph = {a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors()} for a in mol.GetAtoms()}
    if has_spiro_union_shape(mol, set(atoms)):
        best = min(spiro_union_numberings(mol, graph, atoms), key=lambda n: (n.pre_key, n.unsat_key))
        spiro_locants, citation = tuple(best.spiro_key[1]), tuple(best.spiro_key[2])
        position = {a: loc.key[1::-1] for a, loc in best.position_of.items()}
    else:
        best = min(
            _von_baeyer(mol, atoms),
            key=lambda n: (sorted(n.position_of[a] for a in spiro_atoms), n.pre_key, n.unsat_key),
        )
        position = {a: (0, loc) for a, loc in best.position_of.items()}
        spiro_locants = citation = tuple(sorted((loc, 0, "") for a, (_, loc) in position.items() if a in spiro_atoms))

    def first(i):
        own = [position[a] for a in comps[i]["atoms"] if a not in spiro_atoms]
        return min(own or [position[a] for a in comps[i]["atoms"]])

    seniority = [_parent_key(mol, frozenset(c["atoms"]), as_carbon=_replacement_atoms(mol, c)) for c in comps]
    return (
        -sum(len(o) - 1 for o in owners.values() if len(o) > 1),
        0 if saturated else 1,
        spiro_locants if saturated else (),
        tuple(sorted(seniority)),
        tuple(seniority[i] for i in sorted(range(len(comps)), key=first)),
        citation,
    )


_TYPE_KEYS = {
    "phane": _phane_key,
    "spiro": _spiro_key,
    "fused": _fused_key,
    "bridged_fused": _bridged_key,
    "von_baeyer": _von_baeyer_key,
}
