"""Metallacycles (P-69.4, `tmp/bluebook/P6a.txt` lines 8873-8940): a
monocyclic ring with one Group 3-12 or f-block metal (and optional hetero
atoms) replacing ring carbons, named by skeletal replacement
('1-sila-2-ferracyclopentane', '1-iridabenzene'); hydrido, ylidene, ionic
(-ium/-ide), suffix-group and stereo forms. Polycycles: _metallapolycycle.py.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._coordination import _net_charge, collect_ligands
from ._metal_pair import _brackets
from ._numerals import alkane_name, multiplying_prefix
from ._ring_extras import add_ligands, check_charge, ionic_stem, ring_ligand_entries, ring_stereo, stereo_prefix
from ._ring_groups import add_n_entries, is_carboxyl_bond, is_exocyclic_oxo, principal_kind, ring_substituents, with_suffix
from ._substituents import format_substituent_prefixes

_METAL_A_PREFIXES = {
    "Sc": "scanda", "Y": "ytta", "La": "lanthana", "Ce": "cera", "Pr": "praseodyma", "Nd": "neodyma",
    "Pm": "prometha", "Sm": "samara", "Eu": "europa", "Gd": "gadolina", "Tb": "terba", "Dy": "dysprosa",
    "Ho": "holma", "Er": "erba", "Tm": "thula", "Yb": "ytterba", "Lu": "luteta",
    "Ac": "actina", "Th": "thora", "Pa": "protactina", "U": "urana", "Np": "neptuna", "Pu": "plutona",
    "Am": "america", "Cm": "cura",
    "Ti": "titana", "Zr": "zircona", "Hf": "hafna",
    "V": "vanada", "Nb": "niba", "Ta": "tantala",
    "Cr": "chroma", "Mo": "molybda", "W": "tungsta",
    "Mn": "mangana", "Tc": "techneta", "Re": "rhena",
    "Fe": "ferra", "Ru": "ruthena", "Os": "osma",
    "Co": "cobalta", "Rh": "rhoda", "Ir": "irida",
    "Ni": "nickela", "Pd": "pallada", "Pt": "platina",
    "Cu": "cupra", "Ag": "argenta", "Au": "aura",
    "Zn": "zinca", "Cd": "cadma", "Hg": "mercura",
}

_HETERO_A = {8: "oxa", 16: "thia", 7: "aza", 15: "phospha", 33: "arsa", 14: "sila", 32: "germa", 50: "stanna", 82: "plumba", 5: "bora"}
_HETERO_ORDER = [8, 16, 7, 15, 33, 14, 32, 50, 82, 5]


_DONORS = {7, 8, 15, 16, 33}
_DEFAULT_VALENCE = {7: 3, 8: 2, 15: 3, 16: 2, 33: 3}


def _is_sigma(mol, metal_idx, x):
    """A carbon, or a hetero atom bonded by an ordinary covalent bond, as
    opposed to a lone-pair donor."""
    atom = mol.GetAtomWithIdx(x)
    z = atom.GetAtomicNum()
    if z not in _DONORS:
        return True
    bond = mol.GetBondBetweenAtoms(metal_idx, x)
    return (
        bond.GetBondType() != Chem.BondType.DATIVE
        and not atom.GetIsAromatic()
        and atom.GetFormalCharge() == 0
        and atom.GetTotalDegree() <= _DEFAULT_VALENCE[z]
    )


def skeleton_atoms(mol, metal_idx, graph=None, stop=frozenset()):
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
            stack.extend(y for y in graph[x] if y != metal_idx and y not in stop)
        seen |= comp
        sigma = [x for x in graph[metal_idx] if x in comp]
        carbon = any(mol.GetAtomWithIdx(x).GetAtomicNum() == 6 for x in sigma)
        if len(sigma) >= 2 and carbon and all(_is_sigma(mol, metal_idx, x) for x in sigma):
            skeleton |= comp
    return skeleton


def _find_ring_metal(mol):
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_A_PREFIXES and a.IsInRing()]
    if not metals:
        return None
    metal_set = set(metals)
    graph = adjacency(mol)
    skeleton = set(metal_set)
    for m in metals:
        skeleton |= skeleton_atoms(mol, m, graph)
    rings = [r for r in mol.GetRingInfo().AtomRings() if set(r) <= skeleton]
    metal_rings = [r for r in rings if set(r) & metal_set]
    if len(metal_rings) != 1 or len(metal_rings[0]) < 3 or not metal_set <= set(metal_rings[0]):
        return None
    ring = set(metal_rings[0])
    if any(n in metal_set for m in metals for n in graph[m] if n in ring):
        return None
    if len(metals) == 1 and len(ring) == 3:
        a, b = (x for x in ring if x != metals[0])
        if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() >= 2.0:
            return None
    if any(set(r) & ring for r in rings if r is not metal_rings[0]):
        return None
    return tuple(metals), metal_rings[0]


def has_metallacycle_shape(mol) -> bool:
    return _find_ring_metal(mol) is not None


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


def _hetero_text(locants_by_z, metal_locants) -> str:
    """'a' prefixes of the ring: main-group heteroatoms, then the metals in
    the order of Table VI (IR-2.15.3.1 e): `metal_locants` maps symbol -> locants."""
    parts = []
    for z in _HETERO_ORDER:
        locs = sorted(locants_by_z.get(z, []))
        if locs:
            mult = multiplying_prefix(len(locs)) if len(locs) > 1 else ""
            parts.append(f"{','.join(map(str, locs))}-{mult}{_HETERO_A[z]}")
    for symbol in sorted(metal_locants, key=_metal_rank):
        locs = sorted(metal_locants[symbol])
        mult = multiplying_prefix(len(locs)) if len(locs) > 1 else ""
        parts.append(f"{','.join(map(str, locs))}-{mult}{_METAL_A_PREFIXES[symbol]}")
    return "-".join(parts)


def _metal_rank(symbol):
    from ._coordination import _table_vi_rank

    return _table_vi_rank(Chem.GetPeriodicTable().GetAtomicNumber(symbol))


def _ring_double_bonds(mol, graph, ring_set, metal_set):
    doubles = []
    for a, b, *_ in non_single_bonds(mol):
        in_ring = (a in ring_set, b in ring_set)
        if all(in_ring):
            carbons = [mol.GetAtomWithIdx(x).GetAtomicNum() == 6 for x in (a, b)]
            order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
            if order != 2.0 or not any(carbons) or (not all(carbons) and not (metal_set & {a, b})):
                raise UnsupportedStructure("only ring C=C and metal=C double bonds are supported here")
            doubles.append((a, b))
        elif any(in_ring):
            if not is_exocyclic_oxo(mol, a, b, ring_set) and not (metal_set & {a, b}):
                raise UnsupportedStructure("an exocyclic double bond on the ring is not supported here")
        elif not any(m in graph[x] for x in (a, b) for m in metal_set) and not (
            mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(b).GetIsAromatic()
        ) and not is_carboxyl_bond(mol, a, b):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")
    return doubles


def name_metallacycle(mol) -> str:
    metals, ring_atoms = _find_ring_metal(mol)
    metal_set = set(metals)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    if len(metals) == 1:
        charge = check_charge(mol, metals[0], ring_atoms, _net_charge(mol))
    else:
        charge = 0
        if _net_charge(mol) or any(mol.GetAtomWithIdx(i).GetFormalCharge() for i in ring_atoms):
            raise UnsupportedStructure("a charged ring with several metals is not supported here")
    for i in ring_atoms:
        atom = mol.GetAtomWithIdx(i)
        if i not in metal_set and atom.GetAtomicNum() != 6 and atom.GetAtomicNum() not in _HETERO_A:
            raise UnsupportedStructure("this ring heteroatom is not supported here")
        if atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic ring atoms are not supported here")
    graph = adjacency(mol)
    ring_set = set(ring_atoms)
    doubles_bonds = _ring_double_bonds(mol, graph, ring_set, metal_set)

    size = len(ring_atoms)
    ligand_entries = {}
    for m in metals:
        counts, simple_labels, organic, neutral, _ = collect_ligands(mol, mol.GetAtomWithIdx(m), graph, skip=ring_set)
        ligand_entries[m] = ring_ligand_entries(counts, simple_labels, neutral)
    principal = principal_kind(mol, graph, ring_set, metals[0])

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
                if atom not in metal_set and z != 6:
                    by_z.setdefault(z, []).append(locant[atom])
            metal_locants: dict = {}
            for m in metals:
                metal_locants.setdefault(mol.GetAtomWithIdx(m).GetSymbol(), []).append(locant[m])
            hetero_locants = sorted([l for v in by_z.values() for l in v] + [locant[m] for m in metals])
            seniority = [tuple(sorted(by_z.get(z, []))) for z in _HETERO_ORDER] + [
                tuple(sorted(metal_locants[sym])) for sym in sorted(metal_locants, key=_metal_rank)
            ]
            grouped: dict = {}
            suffix_locants, n_entries = [], []
            for atom in order:
                if atom in metal_set:
                    continue
                found, n_found, count = ring_substituents(mol, graph, atom, ring_set, principal)
                suffix_locants += [locant[atom]] * count
                n_entries += n_found
                for name, compound in found:
                    entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                    entry["locants"].append(locant[atom])
            for m in metals:
                add_ligands(grouped, ligand_entries[m], locant[m])
            stereo = ring_stereo(mol, ring_atoms, locant)
            for entry in grouped.values():
                entry["locants"].sort()
            all_locants = sorted(l for e in grouped.values() for l in e["locants"])
            citation = [grouped[k]["locants"] for k in sorted(grouped)]
            key = (hetero_locants, seniority, sorted(suffix_locants), doubles, all_locants, citation, [c for _, c in stereo])
            candidates.append((key, by_z, metal_locants, doubles, grouped, suffix_locants, n_entries, stereo))
    if not candidates:
        raise UnsupportedStructure("this ring's double-bond pattern is not supported here")
    _, by_z, metal_locants, doubles, grouped, suffix_locants, n_entries, stereo = min(candidates, key=lambda c: c[0])

    benzene = size == 6 and len(doubles) == 3
    prefixes = _brackets(format_substituent_prefixes(add_n_entries(grouped, n_entries)))
    hetero = _hetero_text(by_z, metal_locants)
    only = next(iter(metal_locants.values()))[0]
    stem = ionic_stem(with_suffix(_ring_stem(size, doubles, benzene), suffix_locants, principal), only, charge)
    return f"{stereo_prefix(stereo)}{prefixes}{'-' if prefixes else ''}{hetero}{stem}"
