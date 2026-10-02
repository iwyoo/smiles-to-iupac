"""Metallacycles (P-69.4, `tmp/bluebook/P6a.txt` lines 8873-8940): a
monocyclic ring with one Group 4-12 metal replacing a ring carbon, named by
skeletal replacement ('1-titanacyclobutane', '1,1-dichloro-2,3,4,5-
tetramethyl-1-platinacyclopenta-2,4-diene'). Ligands on the metal are cited
as detachable prefixes at locant 1; fused, bridged and metallabenzene rings
are out of scope.
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
        else:
            raise UnsupportedStructure("only halogen, alkyl and phenyl ring substituents are supported here")
    return entries


def _ring_stem(size: int, double_locants) -> str:
    base = alkane_name(size)[:-3]
    if not double_locants:
        return f"cyclo{base}ane"
    count = len(double_locants)
    locs = ",".join(map(str, double_locants))
    ene = {1: "ene", 2: "diene", 3: "triene", 4: "tetraene"}[count]
    return f"cyclo{base}{'a' if count > 1 else ''}-{locs}-{ene}"


def name_metallacycle(mol) -> str:
    metal_idx, ring_atoms = _find_ring_metal(mol)
    metal = mol.GetAtomWithIdx(metal_idx)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetIsotope() != 0 for a in mol.GetAtoms()) or any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in ring_atoms):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).GetIsAromatic()
           for i in ring_atoms if i != metal_idx):
        raise UnsupportedStructure("only an all-carbon ring besides the metal is supported here")
    graph = adjacency(mol)
    ring_set = set(ring_atoms)
    for a, b, *_ in non_single_bonds(mol):
        in_ring = (a in ring_set, b in ring_set)
        if any(in_ring):
            if not all(in_ring) or metal_idx in (a, b):
                raise UnsupportedStructure("only ring C=C double bonds are supported here")
        elif not any(metal_idx in graph[x] for x in (a, b)):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")

    start = ring_atoms.index(metal_idx)
    cycle = [ring_atoms[(start + k) % len(ring_atoms)] for k in range(len(ring_atoms))]
    candidates = []
    for order in (cycle, [cycle[0]] + cycle[:0:-1]):
        locant = {atom: i + 1 for i, atom in enumerate(order)}
        doubles = sorted(
            min(locant[b.GetBeginAtomIdx()], locant[b.GetEndAtomIdx()])
            for b in mol.GetBonds()
            if b.GetBondTypeAsDouble() == 2.0 and b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
        )
        grouped: dict = {}
        for atom in order[1:]:
            for name, compound in _substituent_entries(mol, graph, atom, ring_set):
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(locant[atom])
        for entry in grouped.values():
            entry["locants"].sort()
        all_locants = sorted(l for e in grouped.values() for l in e["locants"])
        candidates.append((doubles, all_locants, grouped))
    doubles, _, grouped = min(candidates, key=lambda c: (c[0], c[1]))

    counts, simple_labels, organic, neutral = collect_ligands(mol, metal, graph, skip=ring_set)
    if "hydrido" in counts:
        raise UnsupportedStructure("a hydrido ligand on the metal is not supported here yet")
    for label, n in counts.items():
        name = _HALO_FOR_LIGAND.get(label, label)
        compound = label in neutral or not (label in simple_labels or _is_plain(label))
        entry = grouped.setdefault(name, {"locants": [], "compound": compound})
        entry["locants"] = sorted(entry["locants"] + [1] * n)

    prefix = _METAL_A_PREFIXES[metal.GetSymbol()]
    prefixes = format_substituent_prefixes(grouped)
    return f"{prefixes}{'-' if prefixes else ''}1-{prefix}{_ring_stem(len(ring_atoms), doubles)}"


def _is_plain(name: str) -> bool:
    return not any(ch.isdigit() for ch in name) and "(" not in name
