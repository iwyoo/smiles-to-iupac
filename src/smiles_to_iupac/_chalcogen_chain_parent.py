"""Chalcogen chains with a hydrogen end named on the mononuclear and polynuclear parent hydrides with an -ol, -thiol,
-selenol, -tellurol or peroxol-family suffix (P-68.4.1.1, P-68.4.2.2, P-68.4.2.3): methyldisulfanol,
methyldioxidanethiol, tellanediselenol, disulfaneperoxol. A run of identical atoms is never split between parent, suffix
and prefix; -OH outranks -OOH (P-41), and a chain without any such suffix is left to the hydride namers."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._hetero_prefixes import EXTENDED_PREFIXES, is_functional_carbon
from ._hydride_chain import _alternating_parent
from ._numerals import alkane_name, multiplying_prefix
from ._substituents import format_mononuclear_prefixes, name_branch

_STEM = {8: "oxidane", 16: "sulfane", 34: "selane", 52: "tellane"}
_SYMBOL = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_REPLACEMENT = {16: "thio", 34: "seleno", 52: "telluro"}
_RANK = {8: 0, 16: 1, 34: 2, 52: 3}
_HYDROXY_WORD = {8: "ol", 16: "thiol", 34: "selenol", 52: "tellurol"}
_PEROXOL_ORDER = [
    (8, 8), (16, 8), (34, 8), (52, 8), (8, 16), (8, 34), (8, 52), (16, 16),
    (34, 16), (52, 16), (16, 34), (16, 52), (34, 34), (52, 34), (34, 52), (52, 52),
]
_ARM_ELEMENTS = {6, 9, 17, 35, 53}


def _peroxol_word(first, second):
    if first == second == 8:
        return "peroxol"
    if first == second:
        return f"di{_REPLACEMENT[first]}peroxol"
    replaced = "".join(sorted(_REPLACEMENT[z] for z in (first, second) if z != 8))
    return f"{_SYMBOL[first]}{_SYMBOL[second]}-{replaced}peroxol"


def _chain(mol, graph):
    """The ordered atoms of the one run of two or more chalcogen atoms, or None."""
    chalcogens = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in _STEM}
    runs, seen = [], set()
    for start in sorted(chalcogens):
        if start in seen:
            continue
        run, stack = {start}, [start]
        while stack:
            for n in graph[stack.pop()]:
                if n in chalcogens and n not in run:
                    run.add(n)
                    stack.append(n)
        seen |= run
        if len(run) > 1:
            runs.append(run)
    if len(runs) != 1:
        return None
    run = runs[0]
    ends = [a for a in run if sum(n in run for n in graph[a]) == 1]
    if len(ends) != 2 or any(sum(n in run for n in graph[a]) > 2 for a in run):
        return None
    order, previous = [ends[0]], None
    while len(order) < len(run):
        onward = [n for n in graph[order[-1]] if n in run and n != previous]
        if len(onward) != 1:
            return None
        previous = order[-1]
        order.append(onward[0])
    return order


def _arm_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n != blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _plain_arm(mol, graph, arm, chain_set):
    """A hydrocarbon or halogenated group, possibly interrupted by single chalcogen atoms, that the chain carries as a
    prefix; any characteristic group of its own would outrank the chain."""
    for i in arm:
        a = mol.GetAtomWithIdx(i)
        z = a.GetAtomicNum()
        if a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons():
            return False
        if z in _STEM:
            if a.GetDegree() != 2 or a.GetTotalNumHs() or any(n in chain_set or mol.GetAtomWithIdx(n).GetAtomicNum() in _STEM for n in graph[i]):
                return False
        elif z not in _ARM_ELEMENTS:
            return False
        elif z == 6 and (is_functional_carbon(mol, i) or any(
            b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtom(a).GetAtomicNum() != 6 for b in a.GetBonds()
        )):
            return False
    return True


def _segments(mol, chain):
    segments = []
    for position, atom in enumerate(chain):
        z = mol.GetAtomWithIdx(atom).GetAtomicNum()
        if segments and segments[-1][2] == z:
            segments[-1][1] = position
        else:
            segments.append([position, position, z])
    return segments


def _alternating_parents(segments):
    """(first, last, terminal element, inner element) of each a[ba]n stretch of the chain (P-21.2.3.1)."""
    found = []
    for first in range(len(segments)):
        for last in range(first + 2, len(segments), 2):
            run = segments[first : last + 1]
            if any(s[1] != s[0] for s in run):
                break
            terminal, inner = run[0][2], run[1][2]
            if any(s[2] != (terminal if k % 2 == 0 else inner) for k, s in enumerate(run)):
                break
            if _RANK[inner] < _RANK[terminal]:
                found.append((run[0][0], run[-1][1], terminal, inner))
    return found


def name_chalcogen_chain_parent(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    ):
        return None
    graph = adjacency(mol)
    chain = _chain(mol, graph)
    if chain is None:
        return None
    chain_set = set(chain)
    atoms = [mol.GetAtomWithIdx(i) for i in chain]
    if any(a.IsInRing() for a in atoms) or any(
        b.GetBondTypeAsDouble() != 1.0 for a in atoms for b in a.GetBonds()
    ):
        return None
    if any(a.GetDegree() != 2 or a.GetTotalNumHs() for a in atoms[1:-1]):
        return None
    ends = (atoms[0], atoms[-1])
    arms = {}
    for end in ends:
        outside = [n for n in graph[end.GetIdx()] if n not in chain_set]
        if len(outside) > 1 or end.GetTotalNumHs() + len(outside) != 1 or end.GetDegree() != 1 + len(outside):
            return None
        if outside:
            if mol.GetAtomWithIdx(outside[0]).GetAtomicNum() != 6:
                return None
            arms[end.GetIdx()] = outside[0]
    hydrogen_ends = [e.GetIdx() for e in ends if e.GetIdx() not in arms]
    length = len(chain)
    elements = [a.GetAtomicNum() for a in atoms]
    if not hydrogen_ends or (length == 2 and (arms or elements[0] == elements[1])):
        return None
    covered = set(chain_set)
    for end_idx, root in arms.items():
        arm = _arm_atoms(graph, root, end_idx)
        if arm & covered or not _plain_arm(mol, graph, arm, chain_set):
            return None
        covered |= arm
    if len(covered) != mol.GetNumAtoms():
        return None

    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    organic = bool(arms)
    segments = _segments(mol, chain)
    alternating = _alternating_parents(segments)
    if any(first == 0 and last == length - 1 for first, last, _, _ in alternating):
        return None
    parents = [(s[0], s[1], None, s[2]) for s in segments]
    if not organic:
        parents += alternating
    token = EXTENDED_PREFIXES.set(True)
    candidates = []
    try:
        for first, last, terminal, inner in parents:
            candidate = _candidate(
                mol, graph, chain, elements, arms, halogens, aromatic, first, last, terminal, inner
            )
            if candidate is not None:
                candidates.append(candidate)
    finally:
        EXTENDED_PREFIXES.reset(token)
    if not candidates:
        return None
    return min(candidates, key=lambda c: c[0])[1]


def _suffix_group(elements, side):
    """(rank, word) of the suffix a side of the chain forms when it is one or two atoms ending in hydrogen, else None."""
    if len(side) == 1:
        z = elements[side[0]]
        return (0, _RANK[z]), _HYDROXY_WORD[z]
    if len(side) == 2:
        pair = (elements[side[0]], elements[side[1]])
        return (1, _PEROXOL_ORDER.index(pair)), _peroxol_word(*pair)
    return None


def _candidate(mol, graph, chain, elements, arms, halogens, aromatic, first, last, terminal, inner):
    n = len(chain)
    parent_size = last - first + 1
    left = list(range(first - 1, -1, -1))
    right = list(range(last + 1, n))
    capable = {}
    for name, side, far in (("left", left, chain[0]), ("right", right, chain[-1])):
        if not side or far in arms:
            continue
        group = _suffix_group(elements, side)
        if group is None:
            continue
        if len(side) == 1 and terminal is None and parent_size == 1 and inner == 8:
            # oxidane itself is no parent for an -ol or -thiol suffix (P-68.4.2.3)
            continue
        capable[name] = group
    chosen = []
    if capable:
        best = min(group[0] for group in capable.values())
        chosen = sorted(name for name, group in capable.items() if group[0] == best)
    if not chosen:
        return None
    suffix_rank, suffix_word = capable[chosen[0]]
    if len(chosen) == 2:
        simple = suffix_word == "peroxol" or suffix_word in _HYDROXY_WORD.values()
        suffix_word = f"di{suffix_word}" if simple else f"bis({suffix_word})"
    entries = []
    try:
        for name, side, parent_atom_index in (("left", left, first), ("right", right, last)):
            if name in chosen:
                continue
            if side:
                entry = name_branch(
                    graph, chain[side[0]], chain[parent_atom_index], halogens, aromatic, mol=mol, unsaturated=True
                )
                if entry[0] == "hydroperoxy":
                    entry = (entry[0], True)
                entries.append(entry)
        for end in (0, n - 1):
            root = arms.get(chain[end])
            inside = (end == 0 and first == 0) or (end == n - 1 and last == n - 1)
            if root is not None and inside:
                entries.append(name_branch(graph, root, chain[end], halogens, aromatic, mol=mol, unsaturated=True))
    except UnsupportedStructure:
        return None
    if terminal is None:
        parent = ("" if parent_size == 1 else multiplying_prefix(parent_size)) + _STEM[inner]
    else:
        parent = _alternating_parent(terminal, inner, (parent_size + 1) // 2)
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    if suffix_word[0] in "aeiou" and parent.endswith("e"):
        body = parent[:-1] + suffix_word
    elif suffix_word[0].isupper():
        body = f"{parent}-{suffix_word}"
    elif not prefix and suffix_word.startswith("di") and suffix_word != "diperoxol" and suffix_word not in {
        "di" + w for w in _HYDROXY_WORD.values()
    }:
        body = f"{parent}({suffix_word})"
    else:
        body = parent + suffix_word
    key = (-len(chosen), suffix_rank, _RANK[inner], -parent_size, prefix + body)
    return key, prefix + body


def name_chalcogen_chain_multiplicative(mol):
    """'1,1'-methylenebis(trisulfane)' (P-68.4.1.2): two identical hydrogen-terminated chains of three or more identical
    chalcogen atoms joined through one end atom each by an unbranched saturated carbon chain."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return None
    graph = adjacency(mol)
    chalcogens = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in _STEM}
    runs, seen = [], set()
    for start in sorted(chalcogens):
        if start in seen:
            continue
        run, stack = {start}, [start]
        while stack:
            for n in graph[stack.pop()]:
                if n in chalcogens and n not in run:
                    run.add(n)
                    stack.append(n)
        seen |= run
        runs.append(run)
    if len(runs) != 2 or len(runs[0]) != len(runs[1]) or len(runs[0]) < 3:
        return None
    if mol.GetNumAtoms() != sum(len(r) for r in runs) + sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms()):
        return None
    element = {mol.GetAtomWithIdx(i).GetAtomicNum() for r in runs for i in r}
    if len(element) != 1 or any(
        mol.GetAtomWithIdx(i).IsInRing() or any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetAtomWithIdx(i).GetBonds())
        for r in runs
        for i in r
    ):
        return None
    attachments = []
    for run in runs:
        ends = [i for i in run if sum(n in run for n in graph[i]) == 1]
        carbons = [(i, n) for i in ends for n in graph[i] if n not in run]
        hydrogen = [i for i in ends if mol.GetAtomWithIdx(i).GetTotalNumHs() == 1]
        if len(ends) != 2 or len(carbons) != 1 or len(hydrogen) != 1 or any(
            sum(n in run for n in graph[i]) != 2 or graph[i] != [n for n in graph[i] if n in run] for i in run if i not in ends
        ):
            return None
        attachments.append(carbons[0][1])
    linker, previous, here = [], None, attachments[0]
    while True:
        atom = mol.GetAtomWithIdx(here)
        if atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() != 2 or atom.IsInRing():
            return None
        linker.append(here)
        onward = [n for n in graph[here] if n != previous and mol.GetAtomWithIdx(n).GetAtomicNum() == 6]
        if not onward:
            break
        previous, here = here, onward[0]
    if here != attachments[1] or len(linker) + sum(len(r) for r in runs) != mol.GetNumAtoms():
        return None
    parent = multiplying_prefix(len(runs[0])) + _STEM[next(iter(element))]
    central = "methylene" if len(linker) == 1 else f"({alkane_name(len(linker))}-1,{len(linker)}-diyl)"
    return f"1,1'-{central}bis({parent})"
