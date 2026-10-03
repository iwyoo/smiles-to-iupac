"""Bicyclic, spiro and polycyclic von Baeyer metallacycles (P-69.4,
`tmp/bluebook/P6a.txt` lines 8905-8940: '6,6-di(eta5-cyclopenta-2,4-dien-1-yl)-
6-titanabicyclo[3.2.0]heptane'): one Group 4-12 metal in a carbon skeleton,
named with the 'a' prefix on the bicyclo/tricyclo/spiro parent. Ligand
components are split off the skeleton; Cp rings arrive as [CH-]/[CH]
fragments. Fused (mancude) parents other than anthracene are out of scope.
"""

from rdkit import Chem

from ._bicyclic import _strip_leaves, bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._coordination import _CP_LABEL, _cp_charge, _net_charge, collect_ligands
from ._metal_pair import _brackets
from ._metallacycle import (
    _METAL_A_PREFIXES,
    apply_repeats,
    ring_ligand_entries,
    skeleton_atoms,
)
from ._numerals import alkane_name
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
    if len(metals) != 1:
        return None
    graph = adjacency(complex_mol)
    skeleton = skeleton_atoms(complex_mol, metals[0], graph)
    nbrs = [n for n in graph[metals[0]] if n in skeleton]
    if len(nbrs) < 2 or (len(nbrs) > 2 and any(b in graph[a] for a in nbrs for b in nbrs)):
        return None
    sub, orig = _submol(complex_mol, _ring_core(graph, skeleton, metals[0]))
    candidates = _candidates(sub)
    if candidates is None:
        return None
    return metals[0], graph, orig, candidates


def has_metallapolycycle_shape(mol) -> bool:
    split = _split(mol)
    return split is not None and _analyse(split[0]) is not None


def name_metallapolycycle(mol) -> str:
    complex_mol, cp_count, cp_charge = _split(mol)
    metal_idx, graph, orig, candidates = _analyse(complex_mol)
    metal = complex_mol.GetAtomWithIdx(metal_idx)
    core_atoms = {orig[i] for i in candidates[0][0]}
    if _net_charge(complex_mol) + cp_charge != 0:
        raise UnsupportedStructure("a charged complex is not supported here")
    if any(a.GetIsotope() != 0 for a in complex_mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    for a in core_atoms - {metal_idx}:
        atom = complex_mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("only a carbon skeleton besides the metal is supported here")
    double_bonds = []
    for a, b, *_ in non_single_bonds(complex_mol):
        in_core = (a in core_atoms, b in core_atoms)
        if all(in_core):
            if metal_idx in (a, b) or complex_mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure("only ring C=C double bonds are supported here")
            double_bonds.append((a, b))
        elif any(in_core):
            if not is_exocyclic_oxo(complex_mol, a, b, core_atoms):
                raise UnsupportedStructure("an exocyclic double bond on the ring is not supported here")
        elif not any(metal_idx in graph[x] for x in (a, b)) and not (
            complex_mol.GetAtomWithIdx(a).GetIsAromatic() and complex_mol.GetAtomWithIdx(b).GetIsAromatic()
        ) and not is_carboxyl_bond(complex_mol, a, b):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")

    counts, simple_labels, organic, neutral, _ = collect_ligands(complex_mol, metal, graph, skip=core_atoms)
    if "hydrido" in counts:
        raise UnsupportedStructure("a hydrido ligand on the metal is not supported here yet")
    ligand_entries, repeats = ring_ligand_entries(counts, simple_labels, neutral)
    if cp_count:
        ligand_entries.append((f"({_CP_LABEL})", False, cp_count))

    principal = principal_kind(complex_mol, graph, core_atoms, metal_idx)
    best = None
    for order, parent, outer_key in candidates:
        full_order = [orig[i] for i in order]
        locant = {atom: i + 1 for i, atom in enumerate(full_order)}
        grouped: dict = {}
        suffix_locants, n_entries = [], []
        for atom in full_order:
            if atom == metal_idx:
                continue
            found, n_found, count = ring_substituents(complex_mol, graph, atom, core_atoms, principal)
            suffix_locants += [locant[atom]] * count
            n_entries += n_found
            for name, compound in found:
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(locant[atom])
        for name, compound, n in ligand_entries:
            entry = grouped.setdefault(name, {"locants": [], "compound": compound})
            entry["locants"] += [locant[metal_idx]] * n
        for entry in grouped.values():
            entry["locants"].sort()
        ene = sorted(tuple(sorted((locant[a], locant[b]))) for a, b in double_bonds)
        key = (
            tuple(outer_key),
            locant[metal_idx],
            sorted(suffix_locants),
            ene,
            sorted(l for e in grouped.values() for l in e["locants"]),
            [grouped[k]["locants"] for k in sorted(grouped)],
        )
        if best is None or key < best[0]:
            best = (key, locant[metal_idx], grouped, ene, parent, len(full_order), suffix_locants, n_entries)
    _, metal_locant, grouped, ene, parent, size, suffix_locants, n_entries = best
    prefixes = apply_repeats(_brackets(format_substituent_prefixes(add_n_entries(grouped, n_entries))), repeats, metal_locant)
    if ene:
        base = alkane_name(size)[:-3]
        locs = ",".join(str(lo) if hi - lo == 1 else f"{lo}({hi})" for lo, hi in ene)
        suffix = {1: "ene", 2: "diene", 3: "triene", 4: "tetraene"}.get(len(ene))
        if suffix is None:
            raise UnsupportedStructure("this unsaturation count is not supported here")
        parent = parent[: parent.rindex("]") + 1] + f"{base}{'a' if len(ene) > 1 else ''}-{locs}-{suffix}"
    stem = f"{metal_locant}-{_METAL_A_PREFIXES[metal.GetSymbol()]}{with_suffix(parent, suffix_locants, principal)}"
    return f"{prefixes}{'-' if prefixes else ''}{stem}"
