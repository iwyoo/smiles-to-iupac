"""Assemblies of identical parent cations joined by one multiplying group (P-73.5.1.1, P-15.3): the cationic centres
come from the same parent hydride, so the name is 'junction-(linker)bis(cation)', for example
'4,4'-(ethane-1,2-diyl)bis(1-methylpyridin-1-ium)'. A ring-system cation may join through any ring atom; an acyclic
cation (phosphanium, ylium) joins through its cationic atom.
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, specified_stereo_elements
from ._multiplicative import _Context, _Unit, _arm_chain, _assemble, _build_tree, _components, _find_center, _fragment_key, _ring_systems
from ._multiplicative_groups import classify
from ._multiplicative_linker import DecompositionRejected, name_component
from ._multiplicative_ring import UnitText

_ACYCLIC_CENTRES = {6, 8, 15, 16, 33, 34, 52}
_CHALCOGENS = {8, 16, 34, 52}
_JUNCTION_IODO = re.compile(r"(?:^|(?<=[-(\[{]))(\d+[a-z]?)-iodo")


def _neutral_copy(mol):
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetFormalCharge(0)
    neutral = editable.GetMol()
    neutral.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(neutral)
    return neutral


def _decompose(mol, centres):
    systems = _ring_systems(mol)
    ring_of = {a: i for i, atoms in enumerate(systems) for a in atoms}
    core_of = {}
    cores = []
    for centre in centres:
        idx = centre.GetIdx()
        if idx in ring_of:
            core = set(systems[ring_of[idx]])
        elif idx in {a for a in core_of}:
            continue
        else:
            if centre.GetAtomicNum() not in _ACYCLIC_CENTRES:
                raise UnsupportedStructure("this acyclic cationic centre is not a multiplied parent cation")
            core = {idx}
        if any(a in core_of for a in core):
            continue
        cores.append(core)
        for a in core:
            core_of[a] = len(cores) - 1
    if len(cores) < 2:
        raise UnsupportedStructure("the cationic centres belong to one parent structure")
    rest = {a.GetIdx() for a in mol.GetAtoms()} - set(core_of)
    owner, seen = {}, set()
    for start in sorted(rest):
        if start in seen:
            continue
        component, stack = {start}, [start]
        while stack:
            for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
                if n.GetIdx() in rest and n.GetIdx() not in component:
                    component.add(n.GetIdx())
                    stack.append(n.GetIdx())
        seen |= component
        touched = {
            core_of[n.GetIdx()]
            for a in component
            for n in mol.GetAtomWithIdx(a).GetNeighbors()
            if n.GetIdx() in core_of
        }
        owner[frozenset(component)] = touched
    linkers = [c for c, t in owner.items() if len(t) >= 2]
    if len(linkers) != 1 or owner[linkers[0]] != set(range(len(cores))):
        raise UnsupportedStructure("the parent cations are not joined by one common linking group")
    if any(
        b.GetBeginAtomIdx() in core_of
        and b.GetEndAtomIdx() in core_of
        and core_of[b.GetBeginAtomIdx()] != core_of[b.GetEndAtomIdx()]
        for b in mol.GetBonds()
    ):
        raise UnsupportedStructure("directly bonded cationic parents are not named multiplicatively")
    linker = set(linkers[0])
    units = []
    for i, core in enumerate(cores):
        atoms = set(core)
        for component, touched in owner.items():
            if touched == {i}:
                atoms |= component
        bridges = [
            (a, n.GetIdx())
            for a in core
            for n in mol.GetAtomWithIdx(a).GetNeighbors()
            if n.GetIdx() in linker
        ]
        order = mol.GetBondBetweenAtoms(*bridges[0]).GetBondTypeAsDouble() if len(bridges) == 1 else 0.0
        ylidene = order == 2.0 and len(core) == 1 and mol.GetAtomWithIdx(next(iter(core))).GetAtomicNum() in _CHALCOGENS
        if order != 1.0 and not ylidene:
            raise UnsupportedStructure("a parent cation must join the linking group by one single bond")
        junction, linker_atom = bridges[0]
        if len(core) == 1 and junction not in {c.GetIdx() for c in centres}:
            raise UnsupportedStructure("an acyclic parent cation joins through its cationic atom")
        if len(core) == 1 and mol.GetAtomWithIdx(junction).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(linker_atom).IsInRing():
            raise UnsupportedStructure("carbon cations on a chain belong to the chain parent hydride, not to a multiplied unit")
        units.append(_Unit(("U", i), frozenset(core), frozenset(atoms), junction, linker_atom, _fragment_key(mol, atoms, junction)))
    if len({u.key for u in units}) != 1:
        raise UnsupportedStructure("the parent cations of an assembly must be identical")
    return units, linker


def _backbone(mol, linker, anchors):
    """The linking group without the substituents hanging off it: leaf atoms are peeled until only the paths between
    the anchors and the rings on them remain."""
    remaining = set(linker)
    peel = [a for a in remaining if a not in anchors and sum(1 for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in remaining) <= 1]
    while peel:
        atom = peel.pop()
        if atom not in remaining:
            continue
        remaining.discard(atom)
        for n in mol.GetAtomWithIdx(atom).GetNeighbors():
            m = n.GetIdx()
            if m in remaining and m not in anchors and sum(1 for k in mol.GetAtomWithIdx(m).GetNeighbors() if k.GetIdx() in remaining) <= 1:
                peel.append(m)
    return remaining


def _unit_molecule(mol, unit, replacement):
    editable = Chem.RWMol(mol)
    marker = None
    order = int(mol.GetBondBetweenAtoms(unit.junction, unit.linker_atom).GetBondTypeAsDouble())
    if replacement is None:
        editable.GetAtomWithIdx(unit.junction).SetNoImplicit(True)
        editable.GetAtomWithIdx(unit.junction).SetNumExplicitHs(editable.GetAtomWithIdx(unit.junction).GetTotalNumHs() + order)
    elif order != 1:
        raise UnsupportedStructure("the junction locant of a ylidene-bonded parent cation cannot be read")
    else:
        marker = editable.AddAtom(Chem.Atom(replacement))
        editable.AddBond(unit.junction, marker, Chem.BondType.SINGLE)
        editable.GetAtomWithIdx(unit.junction).SetNoImplicit(True)
        editable.GetAtomWithIdx(unit.junction).SetNumExplicitHs(editable.GetAtomWithIdx(unit.junction).GetTotalNumHs())
    editable.RemoveBond(unit.junction, unit.linker_atom)
    keep = set(unit.atoms) | ({marker} if marker is not None else set())
    for index in sorted(set(range(editable.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(index)
    fragment = editable.GetMol()
    try:
        Chem.SanitizeMol(fragment)
    except Exception as error:
        raise UnsupportedStructure("the parent cation fragment cannot be named on its own") from error
    return Chem.MolToSmiles(fragment)


def _unit_text(mol, unit):
    from ._polyfunctional import _without_iodo
    from .core import smiles_to_iupac

    bare = smiles_to_iupac(_unit_molecule(mol, unit, None))
    if any(ch.isdigit() for ch in bare):
        probe = smiles_to_iupac(_unit_molecule(mol, unit, 53))
        match = _JUNCTION_IODO.search(probe)
        if match is None:
            raise UnsupportedStructure("the junction locant of the parent cation cannot be read")
        locant, text = match.group(1), _without_iodo(probe, match.group(1))
    else:
        locant, text = None, bare
    if len(unit.ring_atoms) == 1 and mol.GetAtomWithIdx(next(iter(unit.ring_atoms))).GetAtomicNum() == 6:
        substituted = bool(re.match(r"[\d(\[]", text)) or any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in unit.atoms)
    else:
        substituted = len(unit.atoms) > len(unit.ring_atoms)
    has_locants = any(ch.isdigit() for ch in text)
    return UnitText(text, locant, substituted or not has_locants, has_locants)


def name_cation_assembly(mol):
    centres = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() != 1 for a in centres):
        raise UnsupportedStructure("an assembly of parent cations needs one fragment of singly charged centres")
    if specified_stereo_elements(mol):
        raise UnsupportedStructure("stereodescriptors in an assembly of parent cations are not supported yet")
    units, linker = _decompose(mol, centres)
    groups = classify(_neutral_copy(mol))
    if groups is None:
        raise UnsupportedStructure("this linking group holds an atom or group the multiplicative namer cannot place")
    tree = _build_tree(mol)
    if tree is None:
        raise UnsupportedStructure("a ring closure through the linking group is not supported yet")
    systems, node_of, _ = tree
    linker_nodes = list({node_of[a] for a in _backbone(mol, linker, {u.linker_atom for u in units})})
    components = _components(mol, linker_nodes, systems)
    found = _find_center(mol, components, units)
    if found is None:
        raise UnsupportedStructure("the linking group has no single center")
    valid, edges, comp_of = found
    center = valid[0]
    ctx = _Context(groups, None, None, [], set(), {})
    kind, atoms = components[center]
    center_edges = edges[center]
    arm_chain = _arm_chain(mol, center, edges, comp_of, center_edges[0], {u.junction for u in units})
    if arm_chain is None:
        raise UnsupportedStructure("the arms of the linking group differ")
    try:
        central = name_component(mol, kind, atoms, center_edges, ctx)
        arm_parts = []
        for cid, toward_center, toward_unit in arm_chain:
            akind, aatoms = components[cid]
            arm_parts.append(
                name_component(mol, akind, aatoms, [toward_center, toward_unit], ctx, directed=(toward_unit[0], toward_center[0]))
            )
    except DecompositionRejected as error:
        raise UnsupportedStructure("this linking group is not named multiplicatively") from error
    unit = _unit_text(mol, units[0])
    located = unit.junction_locant is not None
    if not located:
        unit.junction_locant = 0
    name = _assemble(len(units), unit, central, arm_parts)
    return name if located else name.split("-", 1)[1]
