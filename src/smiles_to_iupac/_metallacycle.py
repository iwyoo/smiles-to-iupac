"""Metallacycles (P-69.4, `tmp/bluebook/P6a.txt` lines 8873-8940): a
monocyclic ring with one Group 4-12 metal (and optional main-group hetero
atoms) replacing ring carbons, named by skeletal replacement
('1-sila-2-ferracyclopentane', '1-iridabenzene'). The Hantzsch-Widman
alternative ('platinole') names the same structures and is not emitted.
Bicyclic and anthracene-type rings: _metallabicycle.py, _metallaanthracene.py;
tricyclic and other fused metallacycles are out of scope.
"""

import re

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, non_single_bonds, plain_phenyl_substituent_atoms
from ._coordination import collect_ligands
from ._metal_pair import _brackets
from ._numerals import alkane_name, multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_METAL_A_PREFIXES = {
    "Ti": "titana",
    "V": "vanada",
    "Cr": "chroma",
    "Fe": "ferra",
    "Co": "cobalta",
    "Ni": "nickela",
    "Ru": "ruthena",
    "Os": "osma",
    "Ir": "irida",
    "Pt": "platina",
}

_HALO_FOR_LIGAND = {"fluorido": "fluoro", "chlorido": "chloro", "bromido": "bromo", "iodido": "iodo"}
_HETERO_A = {8: "oxa", 16: "thia", 7: "aza", 15: "phospha", 33: "arsa", 14: "sila", 32: "germa", 50: "stanna", 82: "plumba", 5: "bora"}
_HETERO_ORDER = [8, 16, 7, 15, 33, 14, 32, 50, 82, 5]


_DONORS = {7, 8, 15, 16, 33}


def skeleton_atoms(mol, metal_idx, graph=None):
    """The metal plus every atom of a metal-bound component that returns to
    the metal through two or more sigma bonds from non-donor atoms (the ring
    skeleton and its substituents), as opposed to ordinary ligands."""
    graph = graph or adjacency(mol)
    skeleton, seen = {metal_idx}, set()
    for start in graph[metal_idx]:
        if start in seen:
            continue
        comp, stack = set(), [start]
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(y for y in graph[x] if y != metal_idx)
        seen |= comp
        sigma = [x for x in graph[metal_idx] if x in comp]
        if len(sigma) >= 2 and all(mol.GetAtomWithIdx(x).GetAtomicNum() not in _DONORS for x in sigma):
            skeleton |= comp
    return skeleton


def _find_ring_metal(mol):
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_A_PREFIXES and a.IsInRing()]
    if len(metals) != 1:
        return None
    metal_idx = metals[0]
    skeleton = skeleton_atoms(mol, metal_idx)
    rings = [r for r in mol.GetRingInfo().AtomRings() if set(r) <= skeleton]
    metal_rings = [r for r in rings if metal_idx in r]
    if len(metal_rings) != 1 or len(metal_rings[0]) < 3:
        return None
    ring = set(metal_rings[0])
    if any(set(r) & ring for r in rings if r is not metal_rings[0]):
        return None
    return metal_idx, metal_rings[0]


def has_metallacycle_shape(mol) -> bool:
    return _find_ring_metal(mol) is not None


def _substituent_entries(mol, graph, ring_atom, ring_set):
    entries = []
    for n in graph[ring_atom]:
        if n in ring_set:
            continue
        atom = mol.GetAtomWithIdx(n)
        if atom.GetAtomicNum() in HALOGEN_PREFIXES:
            entries.append((HALOGEN_PREFIXES[atom.GetAtomicNum()], False))
        elif atom.GetAtomicNum() == 6:
            if plain_phenyl_substituent_atoms(mol, graph, {n}):
                entries.append(("phenyl", False))
            else:
                entries.append(name_branch(graph, n, ring_atom, {}, mol=mol))
        elif atom.GetAtomicNum() == 8 and atom.GetDegree() == 2 and atom.GetTotalNumHs() == 0:
            carbon = next(x for x in graph[n] if x != ring_atom)
            alkyl = name_branch(graph, carbon, n, {}, mol=mol)[0]
            if not alkyl.endswith("yl") or any(ch.isdigit() or ch in "()" for ch in alkyl):
                raise UnsupportedStructure("only a plain alkoxy ring substituent is supported here")
            entries.append((alkyl[:-2] + "oxy", False))
        else:
            raise UnsupportedStructure("only halogen, alkyl, alkoxy and phenyl ring substituents are supported here")
    return entries


def _ring_stem(size: int, double_locants, benzene: bool) -> str:
    if benzene:
        return "benzene"
    base = alkane_name(size)[:-3]
    if not double_locants:
        return f"cyclo{base}ane"
    count = len(double_locants)
    locs = ",".join(map(str, double_locants))
    ene = {1: "ene", 2: "diene", 3: "triene", 4: "tetraene"}[count]
    return f"cyclo{base}{'a' if count > 1 else ''}-{locs}-{ene}"


def _hetero_text(locants_by_z, metal_locant, metal_prefix) -> str:
    parts = []
    for z in _HETERO_ORDER:
        locs = sorted(locants_by_z.get(z, []))
        if locs:
            mult = multiplying_prefix(len(locs)) if len(locs) > 1 else ""
            parts.append(f"{','.join(map(str, locs))}-{mult}{_HETERO_A[z]}")
    parts.append(f"{metal_locant}-{metal_prefix}")
    return "-".join(parts)


def _ring_double_bonds(mol, graph, ring_set, metal_idx):
    doubles = []
    for a, b, *_ in non_single_bonds(mol):
        in_ring = (a in ring_set, b in ring_set)
        if all(in_ring):
            carbons = [mol.GetAtomWithIdx(x).GetAtomicNum() == 6 for x in (a, b)]
            order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
            if order != 2.0 or not any(carbons) or (not all(carbons) and metal_idx not in (a, b)):
                raise UnsupportedStructure("only ring C=C and metal=C double bonds are supported here")
            doubles.append((a, b))
        elif any(in_ring):
            raise UnsupportedStructure("an exocyclic double bond on the ring is not supported here")
        elif not any(metal_idx in graph[x] for x in (a, b)) and not (
            mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(b).GetIsAromatic()
        ):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")
    return doubles


def name_metallacycle(mol) -> str:
    metal_idx, ring_atoms = _find_ring_metal(mol)
    metal = mol.GetAtomWithIdx(metal_idx)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetIsotope() != 0 for a in mol.GetAtoms()) or any(
        mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in ring_atoms
    ):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    for i in ring_atoms:
        atom = mol.GetAtomWithIdx(i)
        if i != metal_idx and atom.GetAtomicNum() != 6 and atom.GetAtomicNum() not in _HETERO_A:
            raise UnsupportedStructure("this ring heteroatom is not supported here")
        if atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic ring atoms are not supported here")
    graph = adjacency(mol)
    ring_set = set(ring_atoms)
    doubles_bonds = _ring_double_bonds(mol, graph, ring_set, metal_idx)

    size = len(ring_atoms)
    counts, simple_labels, organic, neutral, _ = collect_ligands(mol, metal, graph, skip=ring_set)
    if "hydrido" in counts:
        raise UnsupportedStructure("a hydrido ligand on the metal is not supported here yet")
    ligand_entries, repeats = ring_ligand_entries(counts, simple_labels, neutral)

    candidates = []
    for offset in range(size):
        cycle = [ring_atoms[(offset + k) % size] for k in range(size)]
        for order in (cycle, [cycle[0]] + cycle[:0:-1]):
            locant = {atom: i + 1 for i, atom in enumerate(order)}
            doubles = []
            for a, b in doubles_bonds:
                low, high = sorted((locant[a], locant[b]))
                if high - low != 1:
                    doubles = None
                    break
                doubles.append(low)
            if doubles is None:
                continue
            doubles.sort()
            by_z: dict = {}
            for atom in order:
                z = mol.GetAtomWithIdx(atom).GetAtomicNum()
                if atom != metal_idx and z != 6:
                    by_z.setdefault(z, []).append(locant[atom])
            hetero_locants = sorted([l for v in by_z.values() for l in v] + [locant[metal_idx]])
            seniority = [tuple(sorted(by_z.get(z, []))) for z in _HETERO_ORDER] + [(locant[metal_idx],)]
            grouped: dict = {}
            for atom in order:
                if atom == metal_idx:
                    continue
                for name, compound in _substituent_entries(mol, graph, atom, ring_set):
                    entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                    entry["locants"].append(locant[atom])
            for name, compound, n in ligand_entries:
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"] += [locant[metal_idx]] * n
            for entry in grouped.values():
                entry["locants"].sort()
            all_locants = sorted(l for e in grouped.values() for l in e["locants"])
            citation = [grouped[k]["locants"] for k in sorted(grouped)]
            key = (hetero_locants, seniority, doubles, all_locants, citation)
            candidates.append((key, by_z, locant[metal_idx], doubles, grouped))
    if not candidates:
        raise UnsupportedStructure("this ring's double-bond pattern is not supported here")
    _, by_z, metal_locant, doubles, grouped = min(candidates, key=lambda c: c[0])

    benzene = size == 6 and len(doubles) == 3
    prefixes = apply_repeats(_brackets(format_substituent_prefixes(grouped)), repeats, metal_locant)
    hetero = _hetero_text(by_z, metal_locant, _METAL_A_PREFIXES[metal.GetSymbol()])
    return f"{prefixes}{'-' if prefixes else ''}{hetero}{_ring_stem(size, doubles, benzene)}"


def ring_ligand_entries(counts, simple_labels, neutral):
    """Metal ligands as ring-name prefixes: halo forms, and a chelating
    ligand cited once without its kappa tag. Returns (entries, repeats):
    `repeats` maps a chelate name to its donor count, so the metal locant
    can be repeated per donor ('9,9-[...]') instead of using 'bis'."""
    entries, repeats = [], {}
    for label, n in counts.items():
        compound = label in neutral or not (label in simple_labels or _is_plain(label))
        name = _HALO_FOR_LIGAND.get(label, label)
        match = re.search(r"-\u03ba(\d+)", label)
        if match and n == 1:
            name = label[: match.start()]
            repeats[name] = int(match.group(1))
        elif match:
            raise UnsupportedStructure("several chelating ligands on a ring metal are not supported here")
        entries.append((name, compound, n))
    return entries, repeats


def apply_repeats(prefixes: str, repeats, metal_locant) -> str:
    for name, k in repeats.items():
        shown = f"[{name}]" if "(" in name else name
        needle = f"{metal_locant}-{shown}"
        if needle in prefixes:
            tail = needle[len(str(metal_locant)):]
            prefixes = prefixes.replace(needle, ",".join([str(metal_locant)] * k) + tail, 1)
    return prefixes


def _is_plain(name: str) -> bool:
    return not any(ch.isdigit() for ch in name) and "(" not in name
