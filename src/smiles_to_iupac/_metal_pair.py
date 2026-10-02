"""Directly bonded Group 13-15 metal pairs (P-69.5.3, `tmp/bluebook/P6a.txt`
lines 8971-8985): the senior metal (As > Sb > Bi > Ge > Sn > Pb > Al > Ga >
In > Tl) is the parent hydride, other metals are '-yl' groups
('germylbismuthane'), also through alkyl chains and para-substituted aryls.
The multiplicative 'plumbanetetrayl' form is out of scope.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, non_single_bonds
from ._group13_hydride import GROUP_13_STEMS, GROUP_14_STEMS, GROUP_15_STEMS
from ._numerals import alkyl_name
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_STEMS = {**GROUP_13_STEMS, **GROUP_14_STEMS, **GROUP_15_STEMS}
_SENIORITY = [33, 51, 83, 32, 50, 82, 13, 31, 49, 81]
_MAX_VALENCE = {**{n: 3 for n in GROUP_13_STEMS}, **{n: 4 for n in GROUP_14_STEMS}, **{n: 3 for n in GROUP_15_STEMS}}


def has_metal_pair_shape(mol) -> bool:
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS) >= 2


def _group_stem(atomic_num: int) -> str:
    stem = _STEMS[atomic_num]
    if atomic_num in GROUP_14_STEMS:
        return stem[:-3] + "yl"
    return stem[:-1] + "yl"


def _format(entries) -> str:
    return format_mononuclear_prefixes(entries) if entries else ""


def _brackets(name: str) -> str:
    return name.replace("(\x01", "[").replace("\x02)", "]").replace("\x01", "[").replace("\x02", "]")


def _mark(name: str) -> str:
    return f"\x01{name}\x02" if "(" in name else name


def _is_carboxy(mol, graph, n):
    oxygens = [x for x in graph[n] if mol.GetAtomWithIdx(x).GetAtomicNum() == 8]
    return (
        len(graph[n]) == 3
        and len(oxygens) == 2
        and sorted(mol.GetBondBetweenAtoms(n, o).GetBondTypeAsDouble() for o in oxygens) == [1.0, 2.0]
        and all(mol.GetAtomWithIdx(o).GetDegree() == 1 for o in oxygens)
    )


def _chain_atoms(graph, root, boundary):
    atoms, stack = set(), [root]
    while stack:
        a = stack.pop()
        if a in atoms:
            continue
        atoms.add(a)
        stack.extend(x for x in graph[a] if x != boundary)
    return atoms


def substituted_aryl_name(mol, graph, root, came_from, seen):
    """Benzene ring attached at `root` (locant 1) whose other positions carry
    halogen, alkyl, carboxy or metal-group substituents."""
    info = mol.GetRingInfo()
    ring = next((r for r in info.AtomRings() if root in r), None)
    if (
        ring is None
        or len(ring) != 6
        or any(not mol.GetAtomWithIdx(i).GetIsAromatic() or mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring)
    ):
        raise UnsupportedStructure("only a benzene-ring aryl group is supported here")
    if any(info.NumAtomRings(i) != 1 for i in ring):
        raise UnsupportedStructure("a fused-ring aryl group is not supported here")
    seen.update(ring)
    start = ring.index(root)
    cycle = [ring[(start + k) % 6] for k in range(6)]
    options = []
    for order in (cycle, [cycle[0]] + cycle[:0:-1]):
        grouped: dict = {}
        for position, atom in enumerate(order[1:], start=2):
            for n in graph[atom]:
                if n in ring:
                    continue
                name, compound = _aryl_substituent(mol, graph, n, atom, seen)
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(position)
        options.append(grouped)

    def key(g):
        return sorted(l for e in g.values() for l in e["locants"]), [g[k]["locants"] for k in sorted(g)]

    grouped = min(options, key=key)
    prefixes = format_substituent_prefixes(grouped)
    return _mark(f"{prefixes}phenyl"), bool(prefixes)


def _aryl_substituent(mol, graph, n, ring_atom, seen):
    atom = mol.GetAtomWithIdx(n)
    z = atom.GetAtomicNum()
    if z in HALOGEN_PREFIXES:
        seen.add(n)
        return HALOGEN_PREFIXES[z], False
    if z in _STEMS:
        return _metal_group_name(mol, graph, n, ring_atom, seen), True
    if z == 6 and not atom.GetIsAromatic():
        if _is_carboxy(mol, graph, n):
            seen.update(_chain_atoms(graph, n, ring_atom))
            return "carboxy", False
        seen.update(_chain_atoms(graph, n, ring_atom))
        return name_branch(graph, n, ring_atom, {}, mol=mol)
    raise UnsupportedStructure("this aryl substituent is not supported here")


def _carbon_group(mol, graph, n, idx, seen):
    """A carbon substituent on `idx` starting at `n`: aryl, alkyl, or an
    unbranched chain carrying exactly one metal group."""
    if mol.GetAtomWithIdx(n).GetIsAromatic():
        return substituted_aryl_name(mol, graph, n, idx, seen)
    branch = _chain_atoms(graph, n, idx)
    if not any(mol.GetAtomWithIdx(a).GetAtomicNum() in _STEMS for a in branch):
        seen.update(branch)
        return name_branch(graph, n, idx, {}, mol=mol)
    chain, prev, cur, linked = [], idx, n, {}
    while True:
        chain.append(cur)
        seen.add(cur)
        onward = [x for x in graph[cur] if x != prev]
        metals = [x for x in onward if mol.GetAtomWithIdx(x).GetAtomicNum() in _STEMS]
        carbons = [x for x in onward if mol.GetAtomWithIdx(x).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(x).GetIsAromatic()]
        if len(metals) + len(carbons) != len(onward) or len(metals) > 1 or len(carbons) > 1:
            raise UnsupportedStructure("only an unbranched alkyl chain carrying one metal group is supported here")
        if metals:
            linked[len(chain)] = (metals[0], cur)
        if not carbons:
            break
        prev, cur = cur, carbons[0]
    if len(linked) != 1:
        raise UnsupportedStructure("only an unbranched alkyl chain carrying one metal group is supported here")
    ((position, (metal, carbon)),) = linked.items()
    group = _metal_group_name(mol, graph, metal, carbon, seen)
    return _mark(f"{position}-({group}){alkyl_name(len(chain))}"), True


def _substituents(mol, graph, idx, came_from, seen):
    entries = []
    for n in graph[idx]:
        if n == came_from:
            continue
        z = mol.GetAtomWithIdx(n).GetAtomicNum()
        if z in _STEMS:
            group = _metal_group_name(mol, graph, n, idx, seen)
            entries.append((group, "\x01" in group))
        elif z == 6:
            entries.append(_carbon_group(mol, graph, n, idx, seen))
        else:
            raise UnsupportedStructure("only metal and carbon substituents are supported here")
    seen.add(idx)
    return entries


def _metal_group_name(mol, graph, idx, came_from, seen) -> str:
    atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
    if mol.GetAtomWithIdx(idx).GetDegree() > _MAX_VALENCE[atomic_num]:
        raise UnsupportedStructure("a metal atom exceeds its valence")
    entries = _substituents(mol, graph, idx, came_from, seen)
    return _mark(_format(entries) + _group_stem(atomic_num))


def name_metal_pair(mol) -> str:
    metals = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS]
    parent_num = next(n for n in _SENIORITY if any(a.GetAtomicNum() == n for a in metals))
    parents = [a for a in metals if a.GetAtomicNum() == parent_num]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetFormalCharge() != 0 or a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if len(parents) != 1:
        from ._metal_chain import name_metal_chain

        return name_metal_chain(mol, adjacency(mol), _STEMS, metals, parent_num, _MAX_VALENCE[parent_num])
    (parent,) = parents
    for a, b, *_ in non_single_bonds(mol):
        atoms = (mol.GetAtomWithIdx(a), mol.GetAtomWithIdx(b))
        if not any(x.GetIsAromatic() or x.GetAtomicNum() == 8 for x in atoms):
            raise UnsupportedStructure("an unsaturated substituent is out of scope here")
    if parent.GetDegree() > _MAX_VALENCE[parent_num]:
        raise UnsupportedStructure("a metal atom exceeds its valence")
    graph = adjacency(mol)
    seen: set[int] = set()
    entries = _substituents(mol, graph, parent.GetIdx(), -1, seen)
    if len(seen) != mol.GetNumAtoms():
        raise UnsupportedStructure("this structure contains atoms outside the supported substituent shapes")
    return _brackets(_format(entries) + _STEMS[parent_num])
