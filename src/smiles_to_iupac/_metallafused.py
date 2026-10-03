"""Metallanaphthalenes, -anthracenes and -phenanthrenes (P-69.4,
`tmp/bluebook/P6a.txt` line ~8920: '9,9-[...]-10H-9-platinaanthracene'): a
fused six-membered carbocycle in which one non-fusion position is a Group
4-12 metal ('a' prefix) and one carbon is the indicated-hydrogen sp3 atom.
Standard parent numbering comes from graph templates matched against the
skeleton; hydro prefixes (more than one sp3 carbon) are out of scope.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._coordination import collect_ligands
from ._metal_pair import _brackets
from ._numerals import multiplying_prefix
from ._metallacycle import _METAL_A_PREFIXES, apply_repeats, ring_ligand_entries
from ._ring_groups import add_n_entries, is_carboxyl_bond, is_exocyclic_oxo, principal_kind, ring_substituents, with_suffix
from ._substituents import format_substituent_prefixes


def _template(labels, edges):
    index = {name: i for i, name in enumerate(labels)}
    adj = {name: set() for name in labels}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return labels, adj


_TEMPLATES = {
    "naphthalene": _template(
        ["1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a"],
        [("1", "2"), ("2", "3"), ("3", "4"), ("4", "4a"), ("4a", "5"), ("5", "6"), ("6", "7"), ("7", "8"),
         ("8", "8a"), ("8a", "1"), ("4a", "8a")],
    ),
    "anthracene": _template(
        ["1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"],
        [("1", "2"), ("2", "3"), ("3", "4"), ("4", "4a"), ("4a", "10"), ("10", "10a"), ("10a", "5"), ("5", "6"),
         ("6", "7"), ("7", "8"), ("8", "8a"), ("8a", "9"), ("9", "9a"), ("9a", "1"), ("4a", "9a"), ("10a", "8a")],
    ),
    "indene": _template(
        ["1", "2", "3", "3a", "4", "5", "6", "7", "7a"],
        [("1", "2"), ("2", "3"), ("3", "3a"), ("3a", "4"), ("4", "5"), ("5", "6"), ("6", "7"), ("7", "7a"),
         ("7a", "1"), ("3a", "7a")],
    ),
    "fluorene": _template(
        ["1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "9a"],
        [("1", "2"), ("2", "3"), ("3", "4"), ("4", "4a"), ("4a", "4b"), ("4b", "5"), ("5", "6"), ("6", "7"),
         ("7", "8"), ("8", "8a"), ("8a", "9"), ("9", "9a"), ("9a", "1"), ("4a", "9a"), ("4b", "8a")],
    ),
    "phenanthrene": _template(
        ["1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "10", "10a"],
        [("1", "2"), ("2", "3"), ("3", "4"), ("4", "4a"), ("4a", "4b"), ("4b", "5"), ("5", "6"), ("6", "7"),
         ("7", "8"), ("8", "8a"), ("8a", "9"), ("9", "10"), ("10", "10a"), ("10a", "1"), ("4a", "10a"),
         ("4b", "8a")],
    ),
}


def _locant_key(position: str):
    digits = "".join(ch for ch in position if ch.isdigit())
    return int(digits), position[len(digits):]


def _mappings(system, graph, labels, adj):
    """Every adjacency-preserving bijection system atoms -> template positions."""
    atoms = sorted(system)
    results = []

    def extend(assigned, used):
        if len(assigned) == len(atoms):
            results.append(dict(assigned))
            return
        atom = next(a for a in atoms if a not in assigned and (not assigned or any(n in assigned for n in graph[a])))
        for pos in labels:
            if pos in used:
                continue
            ok = all(
                (assigned[n] in adj[pos]) == (n in graph[atom])
                for n in assigned
            )
            if ok and all(assigned[n] in adj[pos] for n in graph[atom] if n in assigned):
                assigned[atom] = pos
                used.add(pos)
                extend(assigned, used)
                del assigned[atom]
                used.discard(pos)

    extend({}, set())
    return results


def _find(mol):
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_A_PREFIXES]
    if len(metals) != 1:
        return None
    metal = metals[0]
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if len(r) in (5, 6)]
    start = [r for r in rings if metal in r]
    if len(start) != 1:
        return None
    system, frontier = set(start[0]), [start[0]]
    used = [start[0]]
    changed = True
    while changed:
        changed = False
        for r in rings:
            if r in used:
                continue
            if len(r & system) == 2:
                system |= r
                used.append(r)
                changed = True
    if len(used) < 2 or len(system) not in (9, 10, 13, 14):
        return None
    return metal, system


def has_metallafused_shape(mol) -> bool:
    return _find(mol) is not None and _classify(mol) is not None


def _classify(mol):
    found = _find(mol)
    if found is None:
        return None
    metal, system = found
    graph = adjacency(mol)
    sub_graph = {a: {n for n in graph[a] if n in system} for a in system}
    for name, (labels, adj) in _TEMPLATES.items():
        if len(labels) != len(system):
            continue
        maps = [m for m in _mappings(system, sub_graph, labels, adj) if m[metal] in labels and not m[metal].endswith(("a", "b"))]
        if maps:
            return name, maps
    return None


def name_metallafused(mol) -> str:
    metal_idx, system = _find(mol)
    parent, maps = _classify(mol)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    metal = mol.GetAtomWithIdx(metal_idx)
    graph = adjacency(mol)
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in system - {metal_idx}):
        raise UnsupportedStructure("only a carbon skeleton besides the metal is supported here")
    if any(a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")

    sp3 = [
        i for i in system - {metal_idx}
        if not mol.GetAtomWithIdx(i).GetIsAromatic()
        and not any(
            b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtomIdx(i) in system
            for b in mol.GetAtomWithIdx(i).GetBonds()
        )
    ]
    if len(sp3) % 2 != (len(system) - 1) % 2:
        raise UnsupportedStructure("this sp3 ring-carbon pattern is not a hydro derivative of the mancude parent")
    for a, b, *_ in non_single_bonds(mol):
        if (a in system) != (b in system):
            if not any(metal_idx in graph[x] for x in (a, b)) and not is_exocyclic_oxo(mol, a, b, system):
                raise UnsupportedStructure("an exocyclic double bond on the ring is not supported here")
        elif a not in system and not any(metal_idx in graph[x] for x in (a, b)) and not (
            mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(b).GetIsAromatic()
        ) and not is_carboxyl_bond(mol, a, b):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")

    counts, simple_labels, organic, neutral, _ = collect_ligands(mol, metal, graph, skip=system)
    if "hydrido" in counts:
        raise UnsupportedStructure("a hydrido ligand on the metal is not supported here yet")
    ligand_entries, repeats = ring_ligand_entries(counts, simple_labels, neutral)
    principal = principal_kind(mol, graph, system, metal_idx)

    best = None
    for mapping in maps:
        locant = mapping
        metal_locant = _locant_key(locant[metal_idx])
        grouped: dict = {}
        suffix_locants, n_entries = [], []
        for atom in system - {metal_idx}:
            found, n_found, count = ring_substituents(mol, graph, atom, system, principal)
            suffix_locants += [locant[atom]] * count
            n_entries += n_found
            for name, compound in found:
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(int(locant[atom]))
        for name, compound, n in ligand_entries:
            entry = grouped.setdefault(name, {"locants": [], "compound": compound})
            entry["locants"] += [int(locant[metal_idx])] * n
        for entry in grouped.values():
            entry["locants"].sort()
        sp3_locants = sorted((locant[i] for i in sp3), key=_locant_key)
        indicated = sp3_locants[0] if len(sp3) % 2 else None
        hydro = sp3_locants[1:] if len(sp3) % 2 else sp3_locants
        key = (
            metal_locant,
            _locant_key(indicated) if indicated else (0, ""),
            sorted(suffix_locants, key=_locant_key),
            [_locant_key(h) for h in hydro],
            sorted(l for e in grouped.values() for l in e["locants"]),
            [grouped[k]["locants"] for k in sorted(grouped)],
        )
        if best is None or key < best[0]:
            best = (key, int(locant[metal_idx]), indicated, hydro, grouped, suffix_locants, n_entries)
    _, metal_locant, indicated, hydro, grouped, suffix_locants, n_entries = best
    prefixes = apply_repeats(_brackets(format_substituent_prefixes(add_n_entries(grouped, n_entries))), repeats, metal_locant)
    hydrogen = f"{indicated}H-" if indicated else ""
    hydro_text = f"{','.join(hydro)}-{multiplying_prefix(len(hydro)) if len(hydro) > 1 else ''}hydro-" if hydro else ""
    name = with_suffix(parent, suffix_locants, principal)
    stem = f"{hydrogen}{metal_locant}-{_METAL_A_PREFIXES[metal.GetSymbol()]}{name}"
    return f"{prefixes}{'-' if prefixes else ''}{hydro_text}{stem}"
