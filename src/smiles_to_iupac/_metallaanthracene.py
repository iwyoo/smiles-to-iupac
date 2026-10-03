"""10H-9-metallaanthracenes (P-69.4, `tmp/bluebook/P6a.txt` line ~8920:
'9,9-[...]-10H-9-platinaanthracene'): an anthracene skeleton whose C9 is a
Group 4-12 metal and whose C10 is CH2, with benzo-ring halogen/alkyl
substituents and ligands on the metal. Other fused metallacycles are out
of scope.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._coordination import collect_ligands
from ._metallacycle import _HALO_FOR_LIGAND, _METAL_A_PREFIXES, _is_plain, _substituent_entries
from ._substituents import format_substituent_prefixes


def _system(mol):
    info = mol.GetRingInfo()
    rings = [set(r) for r in info.AtomRings()]
    if len(rings) != 3 or any(len(r) != 6 for r in rings):
        return None
    middle = [r for r in rings if all(r & other for other in rings if other is not r)]
    if len(middle) != 1:
        return None
    mid = middle[0]
    outer = [r for r in rings if r is not mid]
    if outer[0] & outer[1]:
        return None
    metals = [i for i in mid if mol.GetAtomWithIdx(i).GetSymbol() in _METAL_A_PREFIXES]
    if len(metals) != 1:
        return None
    return mid, outer, metals[0]


def has_metallaanthracene_shape(mol) -> bool:
    return _system(mol) is not None


def name_metallaanthracene(mol) -> str:
    mid, outer, metal_idx = _system(mol)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    metal = mol.GetAtomWithIdx(metal_idx)
    graph = adjacency(mol)
    fusion = {i for i in mid if any(i in r for r in outer)}
    others = [i for i in mid if i not in fusion and i != metal_idx]
    if len(fusion) != 4 or len(others) != 1:
        raise UnsupportedStructure("this fused ring arrangement is not supported here")
    x = others[0]
    xa = mol.GetAtomWithIdx(x)
    if xa.GetAtomicNum() != 6 or xa.GetIsAromatic() or xa.GetDegree() != 2 or xa.GetTotalNumHs() != 2:
        raise UnsupportedStructure("only a CH2 at position 10 is supported here")
    system = set(mid) | set().union(*outer)
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in system - {metal_idx}):
        raise UnsupportedStructure("only a carbon skeleton besides the metal is supported here")
    for i in system - {metal_idx, x}:
        if not mol.GetAtomWithIdx(i).GetIsAromatic():
            raise UnsupportedStructure("the benzo rings must be aromatic")
    for a, b, *_ in non_single_bonds(mol):
        if (a in system or b in system) and not (
            mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(b).GetIsAromatic()
        ):
            raise UnsupportedStructure("this unsaturation is not supported here")
        if a not in system and b not in system and not any(metal_idx in graph[v] for v in (a, b)):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")

    counts, simple_labels, organic, neutral, _ = collect_ligands(mol, metal, graph, skip=system)
    if "hydrido" in counts:
        raise UnsupportedStructure("a hydrido ligand on the metal is not supported here yet")
    ligand_entries = [
        (_HALO_FOR_LIGAND.get(label, label), label in neutral or not (label in simple_labels or _is_plain(label)), n)
        for label, n in counts.items()
    ]

    metal_nbrs = [n for n in graph[metal_idx] if n in system]
    best = None
    for first in (0, 1):
        a, b = metal_nbrs[first], metal_nbrs[1 - first]
        ring_a = next(r for r in outer if a in r)
        ring_b = next(r for r in outer if b in r)
        f = next(v for v in graph[x] if v in ring_a)
        g = next(v for v in graph[x] if v in ring_b)

        def path(start, end, ring):
            order, prev, cur = [], end, start
            while True:
                nxt = next(v for v in graph[cur] if v in ring and v != prev and v not in (metal_idx, x))
                if nxt == end:
                    return order
                order.append(nxt)
                prev, cur = cur, nxt

        walk_a = path(a, f, ring_a)
        walk_b = path(g, b, ring_b)
        if len(walk_a) != 4 or len(walk_b) != 4:
            raise UnsupportedStructure("this fused ring arrangement is not supported here")
        locant = {metal_idx: 9, x: 10, a: 9.1, f: 4.1, g: 10.1, b: 8.1}
        for i, atom in enumerate(walk_a, start=1):
            locant[atom] = i
        for i, atom in enumerate(walk_b, start=5):
            locant[atom] = i
        grouped: dict = {}
        for atom in system - {metal_idx, x}:
            for name, compound in _substituent_entries(mol, graph, atom, system):
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(locant[atom])
        for name, compound, n in ligand_entries:
            entry = grouped.setdefault(name, {"locants": [], "compound": compound})
            entry["locants"] += [9] * n
        for entry in grouped.values():
            entry["locants"].sort()
        key = (
            sorted(l for e in grouped.values() for l in e["locants"]),
            [grouped[k]["locants"] for k in sorted(grouped)],
        )
        if best is None or key < best[0]:
            best = (key, grouped)
    prefixes = format_substituent_prefixes(best[1])
    return f"{prefixes}{'-' if prefixes else ''}10H-9-{_METAL_A_PREFIXES[metal.GetSymbol()]}anthracene"
