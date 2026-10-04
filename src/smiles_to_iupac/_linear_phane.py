"""Linear phane names (P-26, P-52.2.5.1): four or more unsubstituted benzene rings, two of them terminal, joined
through single C/O/S/NH nodes with at least seven nodes, e.g. 2,4,6-trioxa-1,7(1),3,5(1,3)-tetrabenzenaheptaphane."""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_numbering_general import _HETERO_RANK
from ._numerals import multiplying_prefix, numerical_term

_PREFIX = {8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza"}
_PATTERN = {1: "1,2", 2: "1,3", 3: "1,4"}


def _benzene_rings(mol):
    rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6 or not all(mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring):
            return None
        rings.append(ring)
    return rings


def _nodes(mol):
    rings = _benzene_rings(mol)
    if rings is None or len(rings) < 4:
        return None
    ring_of = {a: i for i, ring in enumerate(rings) for a in ring}
    if any(mol.GetRingInfo().NumAtomRings(a) != 1 for a in ring_of):
        return None
    attachments = {i: [] for i in range(len(rings))}
    chain = []
    for atom in mol.GetAtoms():
        a = atom.GetIdx()
        if a in ring_of:
            outside = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() not in ring_of or ring_of[n.GetIdx()] != ring_of[a]]
            if outside:
                if len(outside) != 1 or atom.GetTotalNumHs() != 0 or outside[0] in ring_of:
                    return None
                attachments[ring_of[a]].append(a)
            elif atom.GetTotalNumHs() != 1:
                return None
            continue
        if atom.GetFormalCharge() or atom.GetIsotope() or atom.GetDegree() != 2:
            return None
        z = atom.GetAtomicNum()
        if z == 6 and atom.GetTotalNumHs() == 2:
            pass
        elif z in (8, 16, 34, 52) and atom.GetTotalNumHs() == 0:
            pass
        elif z == 7 and atom.GetTotalNumHs() == 1:
            pass
        else:
            return None
        if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
            return None
        chain.append(a)
    if any(len(v) not in (1, 2) for v in attachments.values()) or sum(len(v) == 1 for v in attachments.values()) != 2:
        return None
    start = next(i for i, v in attachments.items() if len(v) == 1)
    path = [("ring", start)]
    first = attachments[start][0]
    previous = first
    current = next(n.GetIdx() for n in mol.GetAtomWithIdx(first).GetNeighbors() if n.GetIdx() not in ring_of)
    while True:
        atom = mol.GetAtomWithIdx(current)
        if current in ring_of:
            ring = ring_of[current]
            path.append(("ring", ring))
            others = [x for x in attachments[ring] if x != current]
            if not others:
                break
            exit_atom = others[0]
            path[-1] = ("ring", ring, current, exit_atom)
            nxt = [n.GetIdx() for n in mol.GetAtomWithIdx(exit_atom).GetNeighbors() if n.GetIdx() not in ring_of or ring_of[n.GetIdx()] != ring]
            previous, current = exit_atom, nxt[0]
        else:
            path.append(("atom", current))
            nxt = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous]
            previous, current = current, nxt[0]
    if len(path) < 7 or len(path) != len(rings) + len(chain):
        return None
    return rings, path


def _local_pattern(ring, entry, exit_atom):
    order = list(ring)
    gap = abs(order.index(entry) - order.index(exit_atom))
    return _PATTERN[min(gap, 6 - gap)]


def has_linear_phane_shape(mol) -> bool:
    return len(Chem.GetMolFrags(mol)) == 1 and _nodes(mol) is not None


def name_linear_phane(mol) -> str:
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        raise UnsupportedStructure("stereo on a linear phane is not supported yet")
    found = _nodes(mol)
    if found is None:
        raise UnsupportedStructure("this structure is not a supported linear phane")
    rings, path = found
    best = None
    for nodes in (path, path[::-1]):
        ring_nodes = [(i + 1, n) for i, n in enumerate(nodes) if n[0] == "ring"]
        hetero = [
            (_HETERO_RANK[mol.GetAtomWithIdx(n[1]).GetSymbol()], i + 1)
            for i, n in enumerate(nodes)
            if n[0] == "atom" and mol.GetAtomWithIdx(n[1]).GetAtomicNum() != 6
        ]
        key = ([loc for loc, _ in ring_nodes], sorted(loc for _, loc in hetero))
        if best is None or key < best[0]:
            best = (key, nodes)
    nodes = best[1]
    patterns = {}
    for i, node in enumerate(nodes):
        if node[0] != "ring":
            continue
        pattern = "1" if len(node) == 2 else _local_pattern(rings[node[1]], node[2], node[3])
        patterns.setdefault(pattern, []).append(i + 1)
    groups = ",".join(
        f"{','.join(map(str, locs))}({pattern})" for pattern, locs in sorted(patterns.items(), key=lambda kv: [int(x) for x in kv[0].split(",")])
    )
    hetero_text = {}
    for i, node in enumerate(nodes):
        if node[0] == "atom":
            z = mol.GetAtomWithIdx(node[1]).GetAtomicNum()
            if z != 6:
                hetero_text.setdefault(z, []).append(i + 1)
    prefixes = "-".join(
        f"{','.join(map(str, locs))}-{multiplying_prefix(len(locs)) if len(locs) > 1 else ''}{_PREFIX[z]}"
        for z, locs in sorted(hetero_text.items(), key=lambda kv: _HETERO_RANK[Chem.GetPeriodicTable().GetElementSymbol(kv[0])])
    )
    count = sum(len(v) for v in patterns.values())
    amplificant = f"{multiplying_prefix(count)}benzena"
    parent = f"{numerical_term(len(nodes))}phane"
    return (prefixes + "-" if prefixes else "") + f"{groups}-{amplificant}{parent}"
