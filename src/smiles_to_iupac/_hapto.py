"""Hapto (eta-n) ligands drawn with a bond from the metal to every
coordinated carbon (P-69.2.4/P-69.2.6, `tmp/bluebook/P6a.txt` lines
8711-8775): full-ring arenes and polyenyls ('eta6-benzene',
'eta5-cyclopenta-2,4-dien-1-yl'), allyl, and partially coordinated
monocyclic/bicyclic polyenes with locants ('[(1,2,5,6-eta)-cycloocta-
1,5-diene]'), plus an eta6-phenyl cited inside a larger ligand name
('[2-(eta6-phenyl)ethanamine]', 'triphenyl(eta6-phenyl)borato'). Bridging
(mu) ligands are not handled.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_ETA = "η"
_ENE = {1: "en", 2: "dien", 3: "trien", 4: "tetraen"}


def _fragment(mol, atoms):
    keep = sorted(atoms)
    rw = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        rw.RemoveAtom(idx)
    frag = rw.GetMol()
    for a in frag.GetAtoms():
        a.SetNoImplicit(False)
    Chem.SanitizeMol(frag)
    return frag, {old: new for new, old in enumerate(keep)}


def _name_fragment(frag):
    from .core import smiles_to_iupac

    return smiles_to_iupac(Chem.MolToSmiles(frag))


def _entries(mol, graph, atom, exclude, allowed):
    out = []
    for n in graph[atom]:
        if n in exclude or n not in allowed:
            continue
        z = mol.GetAtomWithIdx(n).GetAtomicNum()
        if z in HALOGEN_PREFIXES:
            out.append((HALOGEN_PREFIXES[z], False))
        elif z == 6:
            out.append(name_branch(graph, n, atom, {}, mol=mol))
        else:
            raise UnsupportedStructure("only halogen and alkyl substituents are supported on a hapto ligand")
    return out


def _polyenyl(mol, graph, order_candidates, allowed):
    """Radical-type name for an odd ring/chain whose first atom carries the
    free valence: lowest substituent locants over the candidate orders."""
    best = None
    for order in order_candidates:
        grouped: dict = {}
        for position, atom in enumerate(order, start=1):
            for name, compound in _entries(mol, graph, atom, set(order), allowed):
                entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                entry["locants"].append(position)
        key = (sorted(l for e in grouped.values() for l in e["locants"]), [grouped[k]["locants"] for k in sorted(grouped)])
        if best is None or key < best[0]:
            best = (key, grouped)
    prefixes = format_substituent_prefixes(best[1])
    return prefixes, best[1]


def _ring_order(graph, ring_atoms, start, direction):
    ring = list(ring_atoms)
    order = [start]
    prev, cur = None, start
    while len(order) < len(ring):
        nbrs = [v for v in graph[cur] if v in ring_atoms and v != prev]
        nxt = nbrs[direction % len(nbrs)] if len(nbrs) > 1 else nbrs[0]
        order.append(nxt)
        prev, cur = cur, nxt
    return order


def hapto_label(mol, metal, donors_in, atoms):
    graph = adjacency(mol)
    donor_idx = {d.GetIdx() for d in donors_in}
    borate = _borate_label(mol, graph, donor_idx, atoms)
    if borate:
        return borate
    if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
        raise UnsupportedStructure("a charged hapto ligand is not supported yet")
    if not any(
        b.GetBeginAtomIdx() in donor_idx
        and b.GetEndAtomIdx() in donor_idx
        and (b.GetIsAromatic() or b.GetBondTypeAsDouble() >= 2.0)
        for b in mol.GetBonds()
    ):
        raise UnsupportedStructure("a metal bonded to several saturated carbons is not a hapto ligand")
    info = mol.GetRingInfo()
    ring = next((set(r) for r in info.AtomRings() if donor_idx <= set(r) and set(r) <= atoms), None)
    n = len(donor_idx)
    if ring is None:
        if n == 3 and not any(mol.GetAtomWithIdx(i).IsInRing() for i in atoms):
            return _allyl(mol, graph, donor_idx, atoms, n)
        if donor_idx == atoms:
            frag, _ = _fragment(mol, atoms)
            return f"{_ETA}{n}-{_name_fragment(frag)}"
        return _partial(mol, atoms, donor_idx, n)
    size = len(ring)
    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if any(a.GetAtomicNum() != 6 for a in ring_atoms):
        raise UnsupportedStructure("only a carbocyclic hapto ligand is supported")
    if n == size:
        if size == 6 and all(a.GetIsAromatic() for a in ring_atoms):
            frag, _ = _fragment(mol, atoms)
            full = _name_fragment(frag)
            return _phenyl_substituent(mol, graph, ring, atoms, full) or f"{_ETA}6-{full}"
        if size % 2 == 1:
            return _odd_ring(mol, graph, ring, size, atoms)
        frag, _ = _fragment(mol, atoms)
        return f"{_ETA}{n}-{_name_fragment(frag)}"
    return _partial(mol, atoms, donor_idx, n)


def _borate_label(mol, graph, donor_idx, atoms):
    """Tetraorganylborate whose one phenyl ring is the eta6 donor
    ('triphenyl(eta6-phenyl)borato', P-69.2.6); None for any other ligand."""
    boron = [i for i in atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 5]
    if len(boron) != 1 or mol.GetAtomWithIdx(boron[0]).GetFormalCharge() != -1:
        return None
    if sum(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms) != 1 or len(graph[boron[0]]) != 4:
        return None
    ring = next((r for r in mol.GetRingInfo().AtomRings() if donor_idx == set(r) and len(r) == 6), None)
    if ring is None or len([n for n in graph[boron[0]] if n in ring]) != 1:
        return None
    others = [n for n in graph[boron[0]] if n not in ring]
    if len(atoms) != 1 + 6 + sum(len(_branch_atoms(graph, n, boron[0])) for n in others):
        return None
    grouped: dict = {}
    for n in others:
        name, compound = name_branch(graph, n, boron[0], {}, aromatic_atoms={a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}, mol=mol)
        entry = grouped.setdefault(name, {"locants": [], "compound": compound})
        entry["locants"].append(1)
    prefixes = format_substituent_prefixes(grouped, omit_locants=True)
    return f"{prefixes}({_ETA}6-phenyl)borato"


def _branch_atoms(graph, root, parent):
    seen = {root}
    stack = [root]
    while stack:
        for v in graph[stack.pop()]:
            if v != parent and v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def _phenyl_substituent(mol, graph, ring, atoms, full):
    """eta6-phenyl cited inside the name of a larger ligand (P-69.2.6);
    None when the ring is the parent of `full`."""
    ring = list(ring)
    attach = [(a, x) for a in ring for x in graph[a] if x in atoms and x not in ring]
    for a, x in attach:
        side = {x}
        stack = [x]
        while stack:
            for v in graph[stack.pop()]:
                if v in atoms and v not in ring and v not in side:
                    side.add(v)
                    stack.append(v)
        allowed = atoms - side
        orders = [_ring_order(graph, set(ring), a, d) for d in (0, 1)]
        try:
            prefixes, _ = _polyenyl(mol, graph, orders, allowed)
        except UnsupportedStructure:
            continue
        plain = f"{prefixes}phenyl"
        hapto = f"{prefixes}-{_ETA}6-phenyl" if prefixes else f"({_ETA}6-phenyl)"
        at = full.find(plain)
        if at < 0:
            continue
        if at and full[at - 1] not in "([{-":
            raise UnsupportedStructure("an eta6-phenyl among several equivalent substituents is not supported yet")
        return full[:at] + hapto + full[at + len(plain):]
    return None


def _odd_ring(mol, graph, ring, size, atoms):
    if atoms != ring and any(mol.GetAtomWithIdx(i).IsInRing() for i in atoms - ring):
        raise UnsupportedStructure("a ring-substituted hapto ring is not supported yet")
    orders = [_ring_order(graph, ring, s, d) for s in ring for d in (0, 1)]
    prefixes, _ = _polyenyl(mol, graph, orders, atoms)
    count = (size - 1) // 2
    locs = ",".join(str(2 * i) for i in range(1, count + 1))
    base = alkane_name(size)[:-3]
    stem = f"cyclo{base}{'a' if count > 1 else ''}-{locs}-{_ENE[count]}-1-yl"
    return f"{_ETA}{size}-{prefixes}{stem}"


def _allyl(mol, graph, donor_idx, atoms, n):
    if n != 3 or any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in donor_idx):
        raise UnsupportedStructure("this hapto ligand shape is not supported yet")
    ends = [i for i in donor_idx if len(set(graph[i]) & donor_idx) == 1]
    middle = [i for i in donor_idx if len(set(graph[i]) & donor_idx) == 2]
    if len(ends) != 2 or len(middle) != 1 or any(mol.GetAtomWithIdx(i).IsInRing() for i in atoms):
        raise UnsupportedStructure("this hapto ligand shape is not supported yet")
    orders = [[ends[0], middle[0], ends[1]], [ends[1], middle[0], ends[0]]]
    prefixes, _ = _polyenyl(mol, graph, orders, atoms)
    if not prefixes and atoms == donor_idx:
        return f"{_ETA}3-allyl"
    return f"{_ETA}3-{prefixes}prop-2-en-1-yl"


def _partial(mol, atoms, donor_idx, n):
    frag, index = _fragment(mol, atoms)
    donors = {index[i] for i in donor_idx}
    ring_info = frag.GetRingInfo()
    name = _name_fragment(frag)
    doubles = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in frag.GetBonds() if b.GetBondTypeAsDouble() == 2.0
    ]
    orders = []
    if ring_info.NumRings() == 1 and frag.GetNumAtoms() == len(ring_info.AtomRings()[0]):
        graph = adjacency(frag)
        ring = set(range(frag.GetNumAtoms()))
        orders = [_ring_order(graph, ring, s, d) for s in ring for d in (0, 1)]
    else:
        from ._bicyclic import find_bicyclic_core, iter_bicyclic_numberings

        core = find_bicyclic_core(frag)
        if core is None or frag.GetNumAtoms() != 2 + sum(len(b) for b in core[2]):
            raise UnsupportedStructure("this partially coordinated ligand is not supported yet")
        orders = list(iter_bicyclic_numberings(core))
    best = None
    for order in orders:
        locant = {atom: i + 1 for i, atom in enumerate(order)}
        ene = sorted(min(locant[a], locant[b]) for a, b in doubles)
        wrapped = sum(abs(locant[a] - locant[b]) > 1 for a, b in doubles)
        key = (wrapped, ene, sorted(locant[d] for d in donors))
        if best is None or key < best:
            best = key
    best = best[1:]
    locs = ",".join(str(x) for x in best[1])
    return f"[({locs}-{_ETA})-{name}]"
