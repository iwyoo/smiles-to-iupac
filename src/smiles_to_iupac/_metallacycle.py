"""Metallacycles (P-69.4, `tmp/bluebook/P6a.txt` lines 8873-8940): a
monocyclic ring with one Group 4-12 metal (and optional main-group hetero
atoms) replacing ring carbons, named by skeletal replacement
('1-sila-2-ferracyclopentane', '1-iridabenzene'). The Hantzsch-Widman
alternative ('platinole') names the same structures and is not emitted.
Fused and bridged rings are out of scope.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, non_single_bonds, plain_phenyl_substituent_atoms
from ._coordination import collect_ligands
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


def _find_ring_metal(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    metals = [idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetSymbol() in _METAL_A_PREFIXES]
    if len(metals) != 1 or len(ring_atoms) < 3:
        return None
    return metals[0], ring_atoms


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
        elif not any(metal_idx in graph[x] for x in (a, b)):
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
    ligand_entries = []
    for label, n in counts.items():
        compound = label in neutral or not (label in simple_labels or _is_plain(label))
        ligand_entries.append((_HALO_FOR_LIGAND.get(label, label), compound, n))

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
    prefixes = format_substituent_prefixes(grouped)
    hetero = _hetero_text(by_z, metal_locant, _METAL_A_PREFIXES[metal.GetSymbol()])
    return f"{prefixes}{'-' if prefixes else ''}{hetero}{_ring_stem(size, doubles, benzene)}"


def _is_plain(name: str) -> bool:
    return not any(ch.isdigit() for ch in name) and "(" not in name
