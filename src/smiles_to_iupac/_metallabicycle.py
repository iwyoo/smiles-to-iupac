"""Saturated bicyclic metallacycles (P-69.4, `tmp/bluebook/P6a.txt` line
8905: '6,6-di(eta5-cyclopenta-2,4-dien-1-yl)-6-titanabicyclo[3.2.0]heptane'):
one Group 4-12 ring metal in a von Baeyer bicyclic carbon skeleton, named
with the 'a' prefix on 'bicyclo[x.y.z]alkane'. Cp rings arrive as separate
[CH-]/[CH] fragments. Unsaturated, polycyclic and hetero-atom rings are
out of scope.
"""

from rdkit import Chem

from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._coordination import _CP_LABEL, _cp_charge, _net_charge, collect_ligands
from ._metallacycle import _HALO_FOR_LIGAND, _METAL_A_PREFIXES, _is_plain, _substituent_entries
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


def has_metallabicycle_shape(mol) -> bool:
    split = _split(mol)
    if split is None:
        return False
    core = find_bicyclic_core(split[0])
    if core is None:
        return False
    bh1, bh2, bridges = core
    atoms = [bh1, bh2] + [a for b in bridges for a in b]
    return sum(split[0].GetAtomWithIdx(a).GetSymbol() in _METAL_A_PREFIXES for a in atoms) == 1


def name_metallabicycle(mol) -> str:
    complex_mol, cp_count, cp_charge = _split(mol)
    core = find_bicyclic_core(complex_mol)
    bh1, bh2, bridges = core
    core_atoms = {bh1, bh2, *(a for b in bridges for a in b)}
    if _net_charge(complex_mol) + cp_charge != 0:
        raise UnsupportedStructure("a charged complex is not supported here")
    if any(a.GetIsotope() != 0 for a in complex_mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    metal_idx = next(a for a in core_atoms if complex_mol.GetAtomWithIdx(a).GetSymbol() in _METAL_A_PREFIXES)
    metal = complex_mol.GetAtomWithIdx(metal_idx)
    for a in core_atoms - {metal_idx}:
        atom = complex_mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("only a carbon skeleton besides the metal is supported here")
    graph = adjacency(complex_mol)
    for a, b, *_ in non_single_bonds(complex_mol):
        if (a in core_atoms or b in core_atoms) and not any(x == metal_idx for x in (a, b)):
            raise UnsupportedStructure("an unsaturated ring is not supported here")
        if not (a in core_atoms and b in core_atoms) and not any(metal_idx in graph[x] for x in (a, b)):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")
        if a in core_atoms and b in core_atoms:
            raise UnsupportedStructure("an unsaturated ring is not supported here")

    counts, simple_labels, organic, neutral, _ = collect_ligands(complex_mol, metal, graph, skip=core_atoms)
    if "hydrido" in counts:
        raise UnsupportedStructure("a hydrido ligand on the metal is not supported here yet")
    ligand_entries = []
    for label, n in counts.items():
        compound = label in neutral or not (label in simple_labels or _is_plain(label))
        ligand_entries.append((_HALO_FOR_LIGAND.get(label, label), compound, n))
    if cp_count:
        ligand_entries.append((f"({_CP_LABEL})", False, cp_count))

    parent = bicyclic_parent_name(core)
    best = None
    for order in iter_bicyclic_numberings(core):
        locant = {atom: i + 1 for i, atom in enumerate(order)}
        grouped: dict = {}
        for atom in order:
            if atom == metal_idx:
                continue
            for name, compound in _substituent_entries(complex_mol, graph, atom, core_atoms):
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(locant[atom])
        for name, compound, n in ligand_entries:
            entry = grouped.setdefault(name, {"locants": [], "compound": compound})
            entry["locants"] += [locant[metal_idx]] * n
        for entry in grouped.values():
            entry["locants"].sort()
        key = (
            locant[metal_idx],
            sorted(l for e in grouped.values() for l in e["locants"]),
            [grouped[k]["locants"] for k in sorted(grouped)],
        )
        if best is None or key < best[0]:
            best = (key, locant[metal_idx], grouped)
    _, metal_locant, grouped = best
    prefixes = format_substituent_prefixes(grouped)
    stem = f"{metal_locant}-{_METAL_A_PREFIXES[metal.GetSymbol()]}{parent}"
    return f"{prefixes}{'-' if prefixes else ''}{stem}"
