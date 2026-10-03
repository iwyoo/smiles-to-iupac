"""Group 13-15 metal chains and multiplicative metal hydrides (P-69.5.3,
`tmp/bluebook/P6a.txt` lines 8971-8985, P-68.2): an unbranched chain of one
metal element is 'hexaethyldistannane'; one junior metal joining k senior
unsubstituted metal hydrides is 'plumbanetetrayltetrakis(stannane)'.
Branched or mixed-element chains are out of scope.
"""

from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._numerals import multiplying_prefix, numerical_term
from ._substituents import format_substituent_prefixes, name_branch


def _carbon_entries(mol, graph, metal_idx, stems, seen):
    entries = []
    for n in graph[metal_idx]:
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


def name_metal_chain(mol, graph, stems, metals, parent_num, max_valence):
    ids = [m.GetIdx() for m in metals if m.GetAtomicNum() == parent_num]
    if any(m.GetAtomicNum() != parent_num for m in metals):
        return _name_multiplicative(mol, graph, stems, metals, parent_num, max_valence)
    adjacent = {i: [n for n in graph[i] if n in ids] for i in ids}
    ends = [i for i in ids if len(adjacent[i]) <= 1]
    if len(ends) != 2 or any(len(v) > 2 for v in adjacent.values()):
        raise UnsupportedStructure("only an unbranched acyclic metal chain is supported here")
    order = [ends[0]]
    while len(order) < len(ids):
        order.append(next(n for n in adjacent[order[-1]] if n not in order))
    seen = set(ids)
    per_atom = []
    for i in order:
        if mol.GetAtomWithIdx(i).GetDegree() > max_valence:
            raise UnsupportedStructure("a metal atom exceeds its valence")
        per_atom.append(_carbon_entries(mol, graph, i, stems, seen))
    if len(seen) != mol.GetNumAtoms():
        raise UnsupportedStructure("this structure contains atoms outside the supported shapes")
    options = []
    for sequence in (per_atom, per_atom[::-1]):
        grouped: dict = {}
        for locant, entries in enumerate(sequence, start=1):
            for name, compound in entries:
                grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
        options.append(grouped)

    def key(g):
        return sorted(l for e in g.values() for l in e["locants"]), [g[k]["locants"] for k in sorted(g)]

    grouped = min(options, key=key)
    return format_substituent_prefixes(grouped) + numerical_term(len(ids)) + stems[parent_num]


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
