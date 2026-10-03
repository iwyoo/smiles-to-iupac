"""Hapto (eta-n) ligands (P-69.2.4/P-69.2.6): full-ring arenes and polyenyls,
allyl, partial rings with locants, an eta6-phenyl cited inside a larger
ligand name, the borato example; acyclic, hetero-ring and fused-ring ligands
come from _hapto_ext.py and bridging (mu) ligands from `bridge_result`.
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._hapto_ext import _ending, chain_label, fused_label, hetero_ring_label, render, ring_label, valence_gaps
from ._numerals import alkane_name
from ._prefix_groups import PrefixNamer
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
        a.SetFormalCharge(0)
        a.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(frag)
    return frag, {old: new for new, old in enumerate(keep)}


def _name_fragment(frag):
    from .core import smiles_to_iupac

    return smiles_to_iupac(Chem.MolToSmiles(frag))


def _entries(mol, graph, atom, exclude, allowed):
    namer = PrefixNamer(mol, graph)
    return [namer.name(n, atom) for n in graph[atom] if n not in exclude and n in allowed]


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


def _ligand_graph(mol, atoms):
    full = adjacency(mol)
    return {a: [n for n in full[a] if n in atoms] for a in atoms}


def hapto_label(mol, metal, donors_in, atoms, sigma=()):
    graph = _ligand_graph(mol, atoms)
    donor_idx = {d.GetIdx() for d in donors_in}
    sigma = frozenset(d.GetIdx() if hasattr(d, 'GetIdx') else d for d in sigma)
    sigma_orders = {i: mol.GetBondBetweenAtoms(metal.GetIdx(), i).GetBondTypeAsDouble() for i in sigma}
    borate = _borate_label(mol, graph, donor_idx, atoms)
    if borate:
        return borate
    if not any(
        b.GetBeginAtomIdx() in donor_idx
        and b.GetEndAtomIdx() in donor_idx
        and (b.GetIsAromatic() or b.GetBondTypeAsDouble() >= 2.0)
        for b in mol.GetBonds()
    ):
        raise UnsupportedStructure("a metal bonded to several saturated carbons is not a hapto ligand")
    info = mol.GetRingInfo()
    rings = [set(r) for r in info.AtomRings() if set(r) <= atoms]
    donor_rings = [r for r in rings if r & donor_idx]
    n = len(donor_idx)
    if sigma:
        neutral = _neutral_with_kappa(mol, donor_idx, atoms, sigma, sigma_orders)
        if neutral:
            return neutral
    if not donor_rings:
        if n == 3 and donor_idx == atoms and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in donor_idx):
            try:
                return _allyl(mol, graph, donor_idx, atoms, n)
            except UnsupportedStructure:
                pass
        if donor_idx == atoms and all(g == 0 for g in valence_gaps(mol, atoms)[0].values()) and not sigma:
            frag, _ = _fragment(mol, atoms)
            return f"{_ETA}{n}-{_name_fragment(frag)}"
        return render(chain_label(mol, graph, donor_idx, atoms, sigma, sigma_orders))
    system = set(donor_rings[0])
    grown = True
    while grown:
        grown = False
        for r in rings:
            if not r <= system and len(r & system) >= 2:
                system |= r
                grown = True
    single = [r for r in donor_rings if donor_idx <= r]
    multiring = len([r for r in rings if r <= system]) > 1
    if multiring:
        try:
            return render(fused_label(mol, graph, system, donor_idx, atoms, sigma))
        except UnsupportedStructure:
            if not single:
                return _partial(mol, atoms, donor_idx, n)
    if not single:
        raise UnsupportedStructure("donors spread over several rings are not supported here")
    ring = single[0]
    size = len(ring)
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring):
        return render(hetero_ring_label(mol, graph, ring, donor_idx, atoms, sigma))
    if n < size or sigma:
        try:
            return render(ring_label(mol, graph, ring, donor_idx, atoms, sigma))
        except UnsupportedStructure:
            return _partial(mol, atoms, donor_idx, n)
    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if size == 6 and all(a.GetIsAromatic() for a in ring_atoms):
        frag, _ = _fragment(mol, atoms)
        full = _name_fragment(frag)
        return _phenyl_substituent(mol, graph, ring, atoms, full) or f"{_ETA}6-{full}"
    if size % 2 == 1:
        return _odd_ring(mol, graph, ring, size, atoms)
    frag, _ = _fragment(mol, atoms)
    return f"{_ETA}{n}-{_name_fragment(frag)}"


def _kappa_text(mol, sigma):
    symbols = sorted(mol.GetAtomWithIdx(i).GetSymbol() for i in sigma)
    seen: dict[str, int] = {}
    cited = []
    for symbol in symbols:
        cited.append(symbol + "'" * seen.get(symbol, 0))
        seen[symbol] = seen.get(symbol, 0) + 1
    count = str(len(symbols)) if len(symbols) > 1 else ""
    return f"-\u03ba{count}{','.join(cited)}"


def _neutral_with_kappa(mol, donor_idx, atoms, sigma, sigma_orders):
    """Neutral molecule bound by a pi group and by lone-pair donors, named
    as a whole ('(E)-eta2-but-2-enal-kO', IR-10.2.5.1 example 19). Only when
    every unsaturated bond is coordinated, so no eta locants are needed."""
    gaps, kekule = valence_gaps(mol, atoms, sigma_orders)
    if any(gaps.values()) or any(o >= 2 and mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i, o in sigma_orders.items()):
        return None
    if any(
        b.GetBondTypeAsDouble() >= 2.0
        and b.GetBeginAtomIdx() in atoms
        and b.GetEndAtomIdx() in atoms
        and not (b.GetBeginAtomIdx() in donor_idx and b.GetEndAtomIdx() in donor_idx)
        and not (b.GetBeginAtomIdx() in sigma or b.GetEndAtomIdx() in sigma)
        for b in kekule.GetBonds()
    ):
        return None
    for symbol in {mol.GetAtomWithIdx(i).GetSymbol() for i in sigma}:
        if sum(1 for i in atoms if mol.GetAtomWithIdx(i).GetSymbol() == symbol) != sum(
            1 for i in sigma if mol.GetAtomWithIdx(i).GetSymbol() == symbol
        ):
            raise UnsupportedStructure("the kappa donor needs a locant in this ligand name")
    frag, _ = _fragment(mol, atoms)
    name = _name_fragment(frag)
    stereo = re.match(r"^\(([\dEZRSrsez,' ]+)\)-", name)
    head = ""
    if stereo:
        head, name = stereo.group(0), name[stereo.end():]
    return f"{head}{_ETA}{len(donor_idx)}-{name}{_kappa_text(mol, sigma)}"


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
            multiplier = next((m for m in ("di", "tri", "tetra") if full[:at].endswith(m)), None)
            start = at - len(multiplier) if multiplier else at
            if prefixes or not multiplier or (start and full[start - 1] not in "([{-"):
                raise UnsupportedStructure("this eta6-phenyl among equivalent substituents is not supported yet")
            rest = {"di": "phenyl", "tri": "diphenyl", "tetra": "triphenyl"}[multiplier]
            return full[:start] + f"({_ETA}6-phenyl){rest}" + full[at + len(plain):]
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
    stem = f"cyclo{base}{'a' if count > 1 else ''}-{locs}-{_ENE[count]}-1-{_ending(mol, ring)}"
    return f"{prefixes}{'-' if prefixes else ''}{_ETA}{size}-{stem}"


def _allyl(mol, graph, donor_idx, atoms, n):
    if n != 3 or any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in donor_idx):
        raise UnsupportedStructure("this hapto ligand shape is not supported yet")
    ends = [i for i in donor_idx if len(set(graph[i]) & donor_idx) == 1]
    middle = [i for i in donor_idx if len(set(graph[i]) & donor_idx) == 2]
    if len(ends) != 2 or len(middle) != 1 or any(mol.GetAtomWithIdx(i).IsInRing() for i in atoms):
        raise UnsupportedStructure("this hapto ligand shape is not supported yet")
    orders = [[ends[0], middle[0], ends[1]], [ends[1], middle[0], ends[0]]]
    prefixes, _ = _polyenyl(mol, graph, orders, atoms)
    ending = _ending(mol, donor_idx)
    if not prefixes and atoms == donor_idx and ending == "yl":
        return f"{_ETA}3-allyl"
    return f"{prefixes}{'-' if prefixes else ''}{_ETA}3-prop-2-en-1-{ending}"


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


def bridge_result(mol, donor_idx, atoms):
    """Locant-bearing name of a ligand bound to several metals (mu-eta)."""
    graph = _ligand_graph(mol, atoms)
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= atoms]
    donor_rings = [r for r in rings if r & donor_idx]
    if not donor_rings:
        return chain_label(mol, graph, donor_idx, atoms)
    system = set(donor_rings[0])
    grown = True
    while grown:
        grown = False
        for r in rings:
            if not r <= system and len(r & system) >= 2:
                system |= r
                grown = True
    if len([r for r in rings if r <= system]) > 1:
        return fused_label(mol, graph, system, donor_idx, atoms)
    ring = donor_rings[0]
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring):
        return hetero_ring_label(mol, graph, ring, donor_idx, atoms)
    return ring_label(mol, graph, ring, donor_idx, atoms)
