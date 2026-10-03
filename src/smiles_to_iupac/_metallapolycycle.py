"""Bicyclic, spiro and polycyclic von Baeyer metallacycles (P-69.4, e.g.
'6,6-di(eta5-cyclopenta-2,4-dien-1-yl)-6-titanabicyclo[3.2.0]heptane'): one
metal in a carbon skeleton named with the 'a' prefix on the bicyclo/tricyclo/
spiro parent, with ene, suffix groups, substituent rings and stereodescriptors.
Ligand components are split off the skeleton; Cp rings arrive as fragments.
"""

from rdkit import Chem

from ._bicyclic import _strip_leaves, bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._coordination import _CP_LABEL, _cp_charge, _net_charge, collect_ligands
from ._metal_pair import _brackets
from ._metallacycle import (
    FREE_VALENCE_PROP,
    _HETERO_A,
    _HETERO_ORDER,
    _METAL_A_PREFIXES,
    _hetero_text,
    _metal_rank,
    skeleton_atoms,
)
from ._numerals import alkane_name
from ._ring_extras import add_ligands, check_charge, ionic_stem, ring_ligand_entries, ring_stereo, stereo_prefix
from ._ring_groups import add_n_entries, is_carboxyl_bond, is_exocyclic_oxo, principal_kind, ring_substituents, with_suffix
from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates
from ._spiro import find_monospiro_atom, iter_monospiro_numberings
from ._substituents import format_substituent_prefixes


def _split(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    complexes = [f for f in frags if any(a.GetSymbol() in _METAL_A_PREFIXES and a.IsInRing() for a in f.GetAtoms())]
    if len(complexes) != 1:
        return None
    others = [f for f in frags if f is not complexes[0]]
    cp_charges = [_cp_charge(f) for f in others]
    if any(c is None for c in cp_charges):
        return None
    return complexes[0], len(others), sum(cp_charges)


def _submol(mol, atoms):
    keep = sorted(atoms)
    rw = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        rw.RemoveAtom(idx)
    sub = rw.GetMol()
    sub.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(sub)
    return sub, keep


def _candidates(sub):
    graph = adjacency(sub)
    core = _strip_leaves(graph)
    if not core:
        return None
    rings = sum(len(v) for v in core.values()) // 2 - len(core) + 1
    if rings == 2:
        bicyclic = find_bicyclic_core(sub)
        if bicyclic is not None:
            parent = bicyclic_parent_name(bicyclic)
            return [(order, parent, ()) for order in iter_bicyclic_numberings(bicyclic)]
        spiro = find_monospiro_atom(sub)
        if spiro is not None:
            return [(order, parent, ()) for parent, order in iter_monospiro_numberings(sub, spiro)]
        return None
    if rings >= 3:
        polycyclic = find_polycyclic_core(sub, rings)
        if polycyclic is not None:
            return list(iter_polycyclic_candidates(polycyclic, rings))
    return None


def _ring_core(graph, skeleton, start):
    """Atoms reachable from `start` over non-bridge bonds: the ring skeleton
    without substituent rings that hang off it through a single bond."""
    order, low, bridges = {}, {}, set()

    def visit(node, parent):
        order[node] = low[node] = len(order)
        for nxt in graph[node]:
            if nxt not in skeleton or nxt == parent:
                continue
            if nxt in order:
                low[node] = min(low[node], order[nxt])
            else:
                visit(nxt, node)
                low[node] = min(low[node], low[nxt])
                if low[nxt] > order[node]:
                    bridges.add(frozenset((node, nxt)))

    visit(start, None)
    core, stack = {start}, [start]
    while stack:
        x = stack.pop()
        for y in graph[x]:
            if y in skeleton and y not in core and frozenset((x, y)) not in bridges:
                core.add(y)
                stack.append(y)
    return core


def _analyse(complex_mol):
    metals = [a.GetIdx() for a in complex_mol.GetAtoms() if a.GetSymbol() in _METAL_A_PREFIXES and a.IsInRing()]
    if not metals:
        return None
    graph = adjacency(complex_mol)
    skeleton = set(metals)
    for m in metals:
        skeleton |= skeleton_atoms(complex_mol, m, graph)
    for m in metals:
        nbrs = [n for n in graph[m] if n in skeleton]
        if len(nbrs) < 2 or (len(nbrs) > 2 and any(b in graph[a] for a in nbrs for b in nbrs)):
            return None
    core = _ring_core(graph, skeleton, metals[0])
    if not set(metals) <= core or any(n in metals for m in metals for n in graph[m] if n in core):
        return None
    sub, orig = _submol(complex_mol, core)
    candidates = _candidates(sub)
    if candidates is None:
        return None
    return tuple(metals), graph, orig, candidates


def has_metallapolycycle_shape(mol) -> bool:
    split = _split(mol)
    return split is not None and _analyse(split[0]) is not None


def name_metallapolycycle(mol) -> str:
    complex_mol, cp_count, cp_charge = _split(mol)
    metals, graph, orig, candidates = _analyse(complex_mol)
    metal_set = set(metals)
    core_atoms = {orig[i] for i in candidates[0][0]}
    charge = _net_charge(complex_mol) + cp_charge
    if len(metals) == 1:
        if charge or complex_mol.GetAtomWithIdx(metals[0]).GetFormalCharge():
            charge = check_charge(complex_mol, metals[0], core_atoms, _net_charge(complex_mol)) if not cp_charge else 0
    elif charge or any(complex_mol.GetAtomWithIdx(i).GetFormalCharge() for i in core_atoms):
        raise UnsupportedStructure("a charged ring with several metals is not supported here")
    if any(a.GetIsotope() != 0 for a in complex_mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    for a in core_atoms - metal_set:
        atom = complex_mol.GetAtomWithIdx(a)
        if (atom.GetAtomicNum() != 6 and atom.GetAtomicNum() not in _HETERO_A) or atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("this ring atom is not supported here")
    double_bonds = []
    for a, b, *_ in non_single_bonds(complex_mol):
        in_core = (a in core_atoms, b in core_atoms)
        if all(in_core):
            if complex_mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure("only ring C=C double bonds are supported here")
            double_bonds.append((a, b))
        elif any(in_core):
            if not is_exocyclic_oxo(complex_mol, a, b, core_atoms) and not (metal_set & {a, b}):
                raise UnsupportedStructure("an exocyclic double bond on the ring is not supported here")
        elif not any(m in graph[x] for x in (a, b) for m in metal_set) and not (
            complex_mol.GetAtomWithIdx(a).GetIsAromatic() and complex_mol.GetAtomWithIdx(b).GetIsAromatic()
        ) and not is_carboxyl_bond(complex_mol, a, b):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")

    ligand_entries = {}
    for m in metals:
        counts, simple_labels, organic, neutral, _ = collect_ligands(complex_mol, complex_mol.GetAtomWithIdx(m), graph, skip=core_atoms)
        ligand_entries[m] = ring_ligand_entries(counts, simple_labels, neutral)
    if cp_count:
        ligand_entries[metals[0]].append((f"({_CP_LABEL})", False, cp_count, 1))

    free = [a.GetIdx() for a in complex_mol.GetAtoms() if a.HasProp(FREE_VALENCE_PROP)]
    principal = "yl" if free else principal_kind(complex_mol, graph, core_atoms, metals[0])
    best = None
    for order, parent, outer_key in candidates:
        full_order = [orig[i] for i in order]
        locant = {atom: i + 1 for i, atom in enumerate(full_order)}
        grouped: dict = {}
        suffix_locants, n_entries = [], []
        by_z: dict = {}
        metal_locants: dict = {}
        for atom in full_order:
            if atom in metal_set:
                metal_locants.setdefault(complex_mol.GetAtomWithIdx(atom).GetSymbol(), []).append(locant[atom])
                continue
            z = complex_mol.GetAtomWithIdx(atom).GetAtomicNum()
            if z != 6:
                by_z.setdefault(z, []).append(locant[atom])
            found, n_found, count = ring_substituents(complex_mol, graph, atom, core_atoms, principal)
            suffix_locants += [locant[atom]] * count
            n_entries += n_found
            for name, compound in found:
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(locant[atom])
        suffix_locants += [locant[i] for i in free]
        for m in metals:
            add_ligands(grouped, ligand_entries[m], locant[m])
        stereo = ring_stereo(complex_mol, core_atoms, locant)
        for entry in grouped.values():
            entry["locants"].sort()
        ene = sorted(tuple(sorted((locant[a], locant[b]))) for a, b in double_bonds)
        hetero_locants = sorted([l for v in by_z.values() for l in v] + [locant[m] for m in metals])
        seniority = [tuple(sorted(by_z.get(z, []))) for z in _HETERO_ORDER] + [
            tuple(sorted(metal_locants[sym])) for sym in sorted(metal_locants, key=_metal_rank)
        ]
        key = (
            tuple(outer_key),
            hetero_locants,
            seniority,
            sorted(suffix_locants),
            ene,
            sorted(l for e in grouped.values() for l in e["locants"]),
            [grouped[k]["locants"] for k in sorted(grouped)],
            [c for _, c in stereo],
        )
        if best is None or key < best[0]:
            best = (key, metal_locants, by_z, grouped, ene, parent, len(full_order), suffix_locants, n_entries, stereo)
    _, metal_locants, by_z, grouped, ene, parent, size, suffix_locants, n_entries, stereo = best
    prefixes = _brackets(format_substituent_prefixes(add_n_entries(grouped, n_entries)))
    if ene:
        base = alkane_name(size)[:-3]
        locs = ",".join(str(lo) if hi - lo == 1 else f"{lo}({hi})" for lo, hi in ene)
        suffix = {1: "ene", 2: "diene", 3: "triene", 4: "tetraene"}.get(len(ene))
        if suffix is None:
            raise UnsupportedStructure("this unsaturation count is not supported here")
        parent = parent[: parent.rindex("]") + 1] + f"{base}{'a' if len(ene) > 1 else ''}-{locs}-{suffix}"
    only = next(iter(metal_locants.values()))[0]
    named = ionic_stem(with_suffix(parent, suffix_locants, principal), only, charge)
    return f"{stereo_prefix(stereo)}{prefixes}{'-' if prefixes else ''}{_hetero_text(by_z, metal_locants)}{named}"
