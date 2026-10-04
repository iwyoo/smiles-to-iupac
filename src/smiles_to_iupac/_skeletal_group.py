"""Skeletal-replacement ('a') substituent groups (P-29.3.2.2, P-29.3.3, P-15.4.3): an unbranched C/N/O/S chain
with four or more heteroatoms, or a saturated heteromonocycle of eleven or more members, entered at a carbon."""

from rdkit import Chem

from ._numerals import alkane_name, multiplying_prefix
from ._common import UnsupportedStructure

_ALLOWED = {6, 7, 8, 16}
_A_PREFIX = {8: "oxa", 16: "thia", 7: "aza"}
_SENIORITY = [8, 16, 7]


def _arm(graph, root, coming_from):
    atoms, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != coming_from and n not in atoms:
                atoms.add(n)
                stack.append(n)
    return atoms


def _locant_text(hetero):
    parts = []
    for z in _SENIORITY:
        locs = sorted(hetero.get(z, []))
        if locs:
            parts.append(f"{','.join(map(str, locs))}-{multiplying_prefix(len(locs)) if len(locs) > 1 else ''}{_A_PREFIX[z]}")
    return "-".join(parts)


def _plain(mol, atoms, coming_from, root):
    if any(
        mol.GetAtomWithIdx(a).GetAtomicNum() not in _ALLOWED or mol.GetAtomWithIdx(a).GetFormalCharge() or mol.GetAtomWithIdx(a).GetIsotope()
        for a in atoms
    ):
        return False
    return all(
        b.GetBondTypeAsDouble() == 1.0 and not b.GetIsAromatic()
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
    ) and mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() == 1.0


def _stereo_free(mol, atoms):
    return all(mol.GetAtomWithIdx(a).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for a in atoms)


def skeletal_chain_group(mol, graph, root, coming_from):
    if mol.GetAtomWithIdx(root).GetAtomicNum() != 6 or mol.GetAtomWithIdx(root).IsInRing():
        return None
    atoms = _arm(graph, root, coming_from)
    if any(mol.GetAtomWithIdx(a).IsInRing() for a in atoms) or not _plain(mol, atoms, coming_from, root) or not _stereo_free(mol, atoms):
        return None
    if any(len([n for n in graph[a] if n in atoms]) > 2 for a in atoms) or len([n for n in graph[root] if n in atoms]) > 1:
        return None
    walk, previous = [root], coming_from
    while True:
        nxt = [n for n in graph[walk[-1]] if n in atoms and n != previous]
        if not nxt:
            break
        previous = walk[-1]
        walk.append(nxt[0])
    zs = [mol.GetAtomWithIdx(a).GetAtomicNum() for a in walk]
    hetero = [i for i, z in enumerate(zs) if z != 6]
    if len(hetero) < 4 or zs[-1] != 6:
        return None
    if any(j - i == 1 for i, j in zip(hetero, hetero[1:])):
        raise UnsupportedStructure("adjacent heteroatoms in a skeletal-replacement group are not supported yet")
    if any(mol.GetAtomWithIdx(a).GetTotalDegree() != len(graph[a]) + mol.GetAtomWithIdx(a).GetTotalNumHs() for a in walk):
        return None
    count = len(walk)
    best = None
    for locate in (lambda i: count - i, lambda i: i + 1):
        found = {z: [locate(i) for i, x in enumerate(zs) if x == z] for z in _SENIORITY}
        key = (sorted(l for v in found.values() for l in v), [sorted(found[z]) for z in _SENIORITY], locate(0))
        if best is None or key < best[0]:
            best = (key, found, locate(0))
    _, found, free = best
    return f"{_locant_text(found)}{alkane_name(count)[:-1]}-{free}-yl", True


def skeletal_ring_group(mol, graph, root, coming_from):
    from ._hetero_chain import has_stereo
    from ._substituents import format_substituent_prefixes, name_branch

    info = mol.GetRingInfo()
    ring = next((list(r) for r in info.AtomRings() if root in r), None)
    if ring is None or info.NumRings() != 1 or len(ring) < 11 or has_stereo(mol):
        return None
    ring_set = set(ring)
    zs = {a: mol.GetAtomWithIdx(a).GetAtomicNum() for a in ring}
    if mol.GetAtomWithIdx(root).GetAtomicNum() != 6 or all(z == 6 for z in zs.values()):
        return None
    if any(z not in _ALLOWED for z in zs.values()) or any(mol.GetAtomWithIdx(a).GetIsAromatic() or mol.GetAtomWithIdx(a).GetFormalCharge() for a in ring):
        return None
    if any(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0 for a in ring for b in graph[a] if b in ring_set):
        return None
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and len(graph[a]) > 2 for a in ring):
        return None
    cycle = [ring[0]]
    while len(cycle) < len(ring):
        cycle.append(next(n for n in graph[cycle[-1]] if n in ring_set and n not in cycle))
    size = len(cycle)
    best = None
    for start in range(size):
        for direction in (1, -1):
            order = [cycle[(start + direction * k) % size] for k in range(size)]
            locant = {a: k + 1 for k, a in enumerate(order)}
            found = {z: [locant[a] for a in order if zs[a] == z] for z in _SENIORITY}
            subs = {}
            for a in order:
                for n in graph[a]:
                    if n not in ring_set and not (a == root and n == coming_from):
                        name, compound = name_branch(graph, n, a, {}, mol=mol)
                        entry = subs.setdefault(name, {"locants": [], "compound": compound})
                        entry["locants"].append(locant[a])
            key = (
                sorted(l for v in found.values() for l in v),
                [found[z] for z in _SENIORITY],
                locant[root],
                sorted(l for e in subs.values() for l in e["locants"]),
                [subs[k]["locants"] for k in sorted(subs)],
            )
            if best is None or key < best[0]:
                best = (key, found, subs, locant[root])
    _, found, subs, free = best
    prefixes = format_substituent_prefixes(subs) if subs else ""
    stem = "cyclo" + alkane_name(size)[:-3] + f"an-{free}-yl"
    return f"{prefixes}{'-' if prefixes else ''}{_locant_text(found)}{stem}", True
