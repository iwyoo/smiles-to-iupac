"""Two identical hydrazine units joined through their N1 atoms by an unbranched carbon chain, named multiplicatively
(P-15.3, P-68.3.1.2.2): '1,1′-ethanediylidenebis(2-phenylhydrazine)'."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_PRIME = "′"
_HALOGENS = (9, 17, 35, 53)


def _carbon_component(graph, mol, start, blocked):
    seen, stack = set(), [start]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(n for n in graph[node] if n not in blocked and mol.GetAtomWithIdx(n).GetAtomicNum() == 6)
    return seen


def _parse(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    nitrogens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    if len(nitrogens) != 4 or any(a.GetFormalCharge() or a.IsInRing() or a.GetIsotope() for a in nitrogens):
        return None
    if any(a.GetAtomicNum() not in (6, 7, *_HALOGENS) for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    pairs = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() == 7 and b.GetAtomicNum() == 7:
            if bond.GetBondTypeAsDouble() != 1.0:
                return None
            pairs.append((a.GetIdx(), b.GetIdx()))
    if len(pairs) != 2 or len({i for pair in pairs for i in pair}) != 4:
        return None
    blocked = {i for pair in pairs for i in pair}
    components = {}
    for index, pair in enumerate(pairs):
        for n in pair:
            for c in graph[n]:
                if c not in blocked:
                    components.setdefault(frozenset(_carbon_component(graph, mol, c, blocked)), []).append((index, n, c))
    linkers = [(comp, hooks) for comp, hooks in components.items() if {u for u, _, _ in hooks} == {0, 1}]
    if len(linkers) != 1 or len(linkers[0][1]) != 2:
        return None
    linker, hooks = linkers[0]
    return graph, pairs, linker, hooks


def _linker_name(mol, graph, linker, hooks):
    """Name of the unbranched, unsubstituted carbon chain whose two ends (or single carbon) are bonded to the units."""
    nitrogens = {h[1] for h in hooks}
    if any(mol.GetAtomWithIdx(a).IsInRing() for a in linker):
        return None
    if any(n not in linker and n not in nitrogens for a in linker for n in graph[a]):
        return None
    internal = {a: [n for n in graph[a] if n in linker] for a in linker}
    if any(len(v) > 2 for v in internal.values()):
        return None
    if any(
        mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble() != 1.0 for a, neighbors in internal.items() for n in neighbors
    ):
        return None
    orders = {mol.GetBondBetweenAtoms(n, c).GetBondTypeAsDouble() for _, n, c in hooks}
    if len(orders) != 1 or next(iter(orders)) not in (1.0, 2.0):
        return None
    double = next(iter(orders)) == 2.0
    length = len(linker)
    if length == 1:
        return "methanediylidene" if double else "methylene"
    ends = {a for a in linker if len(internal[a]) == 1}
    if len(ends) != 2 or {c for _, _, c in hooks} != ends:
        return None
    stem = alkane_name(length)
    if double:
        return f"{stem}diylidene" if length == 2 else f"{stem}-1,{length}-diylidene"
    return f"{stem}-1,{length}-diyl"


def _unit_entries(mol, graph, halogens, pair, hook):
    first, second = hook[1], next(n for n in pair if n != hook[1])
    entries = {}
    for locant, atom, skip in ((1, first, {second, hook[2]}), (2, second, {first})):
        for n in graph[atom]:
            if n not in skip:
                entries.setdefault(locant, []).append(name_branch(graph, n, atom, halogens, mol=mol, unsaturated=True))
    return entries


def _name(mol, found):
    graph, pairs, linker, hooks = found
    linker_text = _linker_name(mol, graph, linker, hooks)
    if linker_text is None:
        return None
    halogens = halogen_substituents(mol)
    units = [_unit_entries(mol, graph, halogens, pairs[i], next(h for h in hooks if h[0] == i)) for i in (0, 1)]
    keys = [sorted((loc, name) for loc, items in unit.items() for name, _ in items) for unit in units]
    if keys[0] != keys[1]:
        return None
    grouped = group_substituents(units[0])
    if not grouped:
        return f"1,1{_PRIME}-{linker_text}dihydrazine"
    return f"1,1{_PRIME}-{linker_text}bis({format_substituent_prefixes(grouped)}hydrazine)"


def has_hydrazine_multiplicative_shape(mol) -> bool:
    found = _parse(mol)
    if found is None:
        return False
    try:
        return _name(mol, found) is not None
    except UnsupportedStructure:
        return False


def name_hydrazine_multiplicative(mol) -> str:
    found = _parse(mol)
    name = _name(mol, found) if found is not None else None
    if name is None:
        raise UnsupportedStructure("this is not two identical hydrazine units joined by a carbon chain")
    return name
