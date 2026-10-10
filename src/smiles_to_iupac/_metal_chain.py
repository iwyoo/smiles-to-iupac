"""Group 13-15 metal chains and multiplicative metal hydrides (P-69.5.3,
the Blue Book, P-68.2): an unbranched chain of one
metal element is 'hexaethyldistannane'; one junior metal joining k senior
unsubstituted metal hydrides is 'plumbanetetrayltetrakis(stannane)'.
Branched or mixed-element chains are out of scope.
"""

from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._numerals import multiplying_prefix, numerical_term
from ._substituents import format_substituent_prefixes, name_branch


def _carbon_entries(mol, graph, metal_idx, stems, seen, only=None):
    entries = []
    for n in graph[metal_idx] if only is None else [only]:
        atom = mol.GetAtomWithIdx(n)
        if atom.GetAtomicNum() in stems:
            continue
        if atom.GetAtomicNum() in HALOGEN_PREFIXES:
            seen.add(n)
            entries.append((HALOGEN_PREFIXES[atom.GetAtomicNum()], False))
        elif atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            from ._metal_pair import _chain_atoms

            branch = _chain_atoms(graph, n, metal_idx)
            if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in branch):
                raise UnsupportedStructure("only plain alkyl substituents are supported on a metal chain")
            seen.update(branch)
            entries.append(name_branch(graph, n, metal_idx, {}, mol=mol))
        else:
            raise UnsupportedStructure("only halogen and alkyl substituents are supported on a metal chain")
    return entries


def _longest_paths(adjacent):
    """Every longest path of the tree of metal atoms (P-44.3: the principal chain has the most skeletal atoms)."""
    best, paths = 0, []
    for start in adjacent:
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            onward = [n for n in adjacent[node] if n not in path]
            if not onward:
                if len(path) > best:
                    best, paths = len(path), []
                if len(path) == best:
                    paths.append(path)
            stack.extend((n, path + [n]) for n in onward)
    return paths


def _entries_on(mol, graph, metal_idx, chain, stems, seen):
    entries = []
    for n in graph[metal_idx]:
        atom = mol.GetAtomWithIdx(n)
        if n in chain:
            continue
        if atom.GetAtomicNum() in stems:
            from ._metal_pair import _chain_atoms

            branch = _chain_atoms(graph, n, metal_idx)
            if branch & chain:
                raise UnsupportedStructure("a metal skeleton that closes a ring is not supported")
            seen.update(branch)
            entries.append(name_branch(graph, n, metal_idx, {}, mol=mol))
        else:
            entries.extend(_carbon_entries(mol, graph, metal_idx, set(), seen, only=n))
    return entries


def name_metal_chain(mol, graph, stems, metals, parent_num, max_valence):
    ids = [m.GetIdx() for m in metals if m.GetAtomicNum() == parent_num]
    if any(m.GetAtomicNum() != parent_num for m in metals):
        return _name_multiplicative(mol, graph, stems, metals, parent_num, max_valence)
    adjacent = {i: [n for n in graph[i] if n in ids] for i in ids}
    if len(ids) > 1 and not all(adjacent.values()):
        raise UnsupportedStructure("only a connected metal skeleton is supported here")
    if sum(len(v) for v in adjacent.values()) // 2 != len(ids) - 1:
        raise UnsupportedStructure("a metal skeleton that closes a ring is not supported")
    if any(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0 for a in ids for b in adjacent[a]):
        raise UnsupportedStructure("a non-single bond within a metal chain is not supported yet (P-68.2)")
    if any(mol.GetAtomWithIdx(i).GetDegree() > max_valence for i in ids):
        raise UnsupportedStructure("a metal atom exceeds its valence")
    options = []
    for path in _longest_paths(adjacent):
        seen = set(ids)
        per_atom = [_entries_on(mol, graph, i, set(path), stems, seen) for i in path]
        if len(seen) != mol.GetNumAtoms():
            raise UnsupportedStructure("this structure contains atoms outside the supported shapes")
        for sequence in (per_atom, per_atom[::-1]):
            grouped: dict = {}
            for locant, entries in enumerate(sequence, start=1):
                for name, compound in entries:
                    grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
            options.append((grouped, len(path)))

    def key(g):
        return sorted(l for e in g.values() for l in e["locants"]), [g[k]["locants"] for k in sorted(g)]

    grouped, length = min(options, key=lambda option: key(option[0]))
    return format_substituent_prefixes(grouped) + numerical_term(length) + stems[parent_num]


def _name_multiplicative(mol, graph, stems, metals, parent_num, max_valence):
    juniors = [m for m in metals if m.GetAtomicNum() != parent_num]
    if len(juniors) != 1:
        raise UnsupportedStructure("only a single central metal is supported in a multiplicative name")
    (center,) = juniors
    seniors = [m for m in metals if m.GetAtomicNum() == parent_num]
    center_nbrs = set(graph[center.GetIdx()])
    if len(seniors) < 2 or {m.GetIdx() for m in seniors} != center_nbrs or mol.GetNumAtoms() != len(seniors) + 1:
        raise UnsupportedStructure("only a central metal joining unsubstituted senior metal hydrides is supported")
    k = len(seniors)
    stem = stems[center.GetAtomicNum()]
    suffix = {2: "diyl", 3: "triyl", 4: "tetrayl"}.get(k)
    if suffix is None:
        raise UnsupportedStructure("this multiplicity is not supported")
    return f"{stem}{suffix}{multiplying_prefix(k, compound=True)}({stems[parent_num]})"
