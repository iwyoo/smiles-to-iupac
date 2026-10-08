"""Acyclic C/N/O/S skeletons of polyethers and polyamines (P-15.4.3, P-51.4,
P-63.2.4): an unbranched chain with four or more separated heteroatoms gets
a skeletal-replacement name ('2,5,8,11-tetraoxadodecane'); otherwise every
heteroatom group off the chosen carbon backbone is collapsed into a named
prefix ('1,2-dimethoxyethane', 'N1-(2-aminoethyl)ethane-1,2-diamine').
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, nonstandard_bonding
from ._numerals import alkane_name, multiplying_prefix
from ._phosphanyl_group import PREFIX_PROP
from ._prefix_groups import PrefixNamer
from ._substituents import format_substituent_prefixes, name_branch

_ALLOWED = {6, 7, 8, 16}
_A_PREFIX = {8: "oxa", 16: "thia", 7: "aza"}
_SENIORITY = [8, 16, 7]


def has_stereo(mol) -> bool:
    return any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )


def _acyclic_single(mol) -> bool:
    return (
        not has_stereo(mol)
        and mol.GetRingInfo().NumRings() == 0
        and len(Chem.GetMolFrags(mol)) == 1
        and all(a.GetAtomicNum() in _ALLOWED and a.GetFormalCharge() == 0 and a.GetIsotope() == 0 for a in mol.GetAtoms())
        and all(b.GetBondTypeAsDouble() == 1.0 for b in mol.GetBonds())
    )


def name_skeletal_chain(mol):
    """Skeletal-replacement name of an unbranched chain with >= 4 heteroatoms, or None."""
    if not _acyclic_single(mol):
        return None
    graph = adjacency(mol)
    if any(len(v) > 2 for v in graph.values()):
        return None
    ends = [a for a, v in graph.items() if len(v) == 1]
    if len(ends) != 2 or any(mol.GetAtomWithIdx(e).GetAtomicNum() != 6 for e in ends):
        return None
    chain, prev = [ends[0]], None
    while len(chain) < mol.GetNumAtoms():
        nxt = next(v for v in graph[chain[-1]] if v != prev)
        prev = chain[-1]
        chain.append(nxt)
    hetero = [i for i, a in enumerate(chain) if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(hetero) < 4:
        return None
    if any(b - a == 1 for a, b in zip(hetero, hetero[1:])):
        raise UnsupportedStructure("adjacent heteroatoms in a chain are not supported here")
    best = None
    for direction in (chain, chain[::-1]):
        zs = [mol.GetAtomWithIdx(a).GetAtomicNum() for a in direction]
        lam = {i + 1: n for i, a in enumerate(direction) if (n := nonstandard_bonding(mol.GetAtomWithIdx(a)))}
        locants = {z: [i + 1 for i, x in enumerate(zs) if x == z] for z in _SENIORITY}
        all_locants = sorted(l for v in locants.values() for l in v)
        key = (all_locants, [locants[z] for z in _SENIORITY], sorted(lam), [-lam[p] for p in sorted(lam)])
        if best is None or key < best[0]:
            best = (key, locants, lam, len(zs))
    _, locants, lam, length = best
    parts = []
    for z in _SENIORITY:
        locs = locants[z]
        if locs:
            mult = multiplying_prefix(len(locs)) if len(locs) > 1 else ""
            parts.append(f"{','.join(_cite(p, lam) for p in locs)}-{mult}{_A_PREFIX[z]}")
    return "-".join(parts) + alkane_name(length)


def _cite(locant, lam):
    return f"{locant}λ{lam[locant]}" if locant in lam else str(locant)


def _components(mol, graph):
    seen, comps = set(), []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        if a.GetAtomicNum() != 6 or i in seen:
            continue
        comp, stack = set(), [i]
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(y for y in graph[x] if mol.GetAtomWithIdx(y).GetAtomicNum() == 6)
        seen |= comp
        comps.append(comp)
    return comps


def _beyond(graph, start, parent):
    seen, stack = {start}, [start]
    while stack:
        for v in graph[stack.pop()]:
            if v != parent and v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


class _Contractor:
    def __init__(self, mol, graph, comps):
        self.mol, self.graph, self.comps = mol, graph, comps
        self.namer = PrefixNamer(mol, graph)
        self.names: dict[int, str] = {}

    def _component_of(self, carbon):
        return next(c for c in self.comps if carbon in c)

    def group(self, x, parent):
        if x not in self.names:
            self.names[x] = self.namer.placeholder(x, parent)
        return self.names[x]


def contract_hetero_groups(mol):
    """`mol` with every heteroatom group off the chosen backbone collapsed
    into a named placeholder prefix, or None when nothing applies."""
    if not _acyclic_single(mol):
        return None
    graph = adjacency(mol)
    hydroxyls = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8 and a.GetDegree() == 1 and a.GetTotalNumHs() == 1
    ]
    if any(a.GetAtomicNum() == 16 and a.GetTotalNumHs() for a in mol.GetAtoms()):
        raise UnsupportedStructure("a thiol among other heteroatom groups is not supported here")
    principal = set(hydroxyls) or {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7}
    comps = _components(mol, graph)
    if not comps:
        return None

    def score(comp):
        attached = sum(1 for c in comp for x in graph[c] if x in principal)
        hetero = sum(1 for c in comp for x in graph[c] if mol.GetAtomWithIdx(x).GetAtomicNum() != 6)
        return (attached, len(comp), hetero)

    backbone = max(comps, key=score)
    contractor = _Contractor(mol, graph, comps)
    outermost = []

    def collect(comp, skip, stay):
        for c in comp:
            for y in graph[c]:
                if mol.GetAtomWithIdx(y).GetAtomicNum() == 6 or y == skip:
                    continue
                if y in stay:
                    for q in graph[y]:
                        if q != c and q not in comp:
                            if mol.GetAtomWithIdx(q).GetAtomicNum() != 6:
                                raise UnsupportedStructure("adjacent heteroatoms are not supported here")
                            collect(contractor._component_of(q), y, set())
                else:
                    outermost.append((y, c))

    collect(backbone, None, principal)
    if not outermost:
        return None
    for y, c in outermost:
        contractor.group(y, c)
    rw = Chem.RWMol(mol)
    drop: set[int] = set()
    for y, c in outermost:
        drop |= _beyond(graph, y, c)
        placeholder = rw.AddAtom(Chem.Atom(53))
        rw.AddBond(c, placeholder, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(placeholder).SetProp(PREFIX_PROP, contractor.names[y])
    for idx in sorted(drop, reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


_HW_STEM = {9: "onane", 10: "ecane"}


def name_hetero_macrocycle(mol):
    """Saturated monocycle of >= 9 members with >= 2 O/S/N heteroatoms and
    only alkyl substituents: '1,4,7-triazonane' (Hantzsch-Widman, 9-10
    members) or '1,4,7,10,13,16-hexaoxacyclooctadecane' ('a' prefixes on a
    cycloalkane, P-22.2.3)."""
    info = mol.GetRingInfo()
    if info.NumRings() != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    ring = list(info.AtomRings()[0])
    size = len(ring)
    if size < 9 or any(a.GetAtomicNum() not in _ALLOWED for a in mol.GetAtoms()):
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetIsAromatic() for a in mol.GetAtoms()) or has_stereo(mol):
        return None
    if any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()):
        return None
    ring_set = set(ring)
    zs = {i: mol.GetAtomWithIdx(i).GetAtomicNum() for i in ring}
    if sum(1 for z in zs.values() if z != 6) < 2:
        return None
    graph = adjacency(mol)
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in range(mol.GetNumAtoms()) if i not in ring_set):
        raise UnsupportedStructure("only alkyl substituents are supported on a heteromacrocycle")
    cycle = [ring[0]]
    while len(cycle) < size:
        cycle.append(next(n for n in graph[cycle[-1]] if n in ring_set and n not in cycle))
    best = None
    for start in range(size):
        for direction in (1, -1):
            order = [cycle[(start + direction * k) % size] for k in range(size)]
            locant = {a: k + 1 for k, a in enumerate(order)}
            hetero = {z: [locant[a] for a in order if zs[a] == z] for z in _SENIORITY}
            lam = {locant[a]: n for a in order if (n := nonstandard_bonding(mol.GetAtomWithIdx(a)))}
            substituents: dict = {}
            for a in order:
                for n in graph[a]:
                    if n not in ring_set:
                        name, compound = name_branch(graph, n, a, {}, mol=mol)
                        entry = substituents.setdefault(name, {"locants": [], "compound": compound})
                        entry["locants"].append(locant[a])
            key = (
                sorted(l for v in hetero.values() for l in v),
                [hetero[z] for z in _SENIORITY],
                sorted(lam),
                [-lam[p] for p in sorted(lam)],
                sorted(l for e in substituents.values() for l in e["locants"]),
                [substituents[k]["locants"] for k in sorted(substituents)],
            )
            if best is None or key < best[0]:
                best = (key, hetero, substituents, lam)
    _, hetero, substituents, lam = best
    prefixes = format_substituent_prefixes(substituents)
    ordered = [z for z in _SENIORITY if hetero[z]]
    all_replaced = len(ordered) == 1 and len(hetero[ordered[0]]) == size
    if size in _HW_STEM:
        a_text = ""
        for z in ordered:
            mult = multiplying_prefix(len(hetero[z])) if len(hetero[z]) > 1 else ""
            if mult.endswith("a") and _A_PREFIX[z][0] in "aeiou":
                mult = mult[:-1]
            a_text += mult + _A_PREFIX[z]
        all_locs = "" if all_replaced and not lam else ",".join(_cite(l, lam) for l in sorted(sum(hetero.values(), []))) + "-"
        name = f"{all_locs}{a_text[:-1]}{_HW_STEM[size]}"
    elif all_replaced and not lam:
        z = ordered[0]
        name = multiplying_prefix(size) + _A_PREFIX[z] + "cyclo" + alkane_name(size)[:-3] + "ane"
    else:
        parts = [
            f"{','.join(_cite(p, lam) for p in sorted(hetero[z]))}-{multiplying_prefix(len(hetero[z])) if len(hetero[z]) > 1 else ''}{_A_PREFIX[z]}"
            for z in ordered
        ]
        name = "-".join(parts) + "cyclo" + alkane_name(size)[:-3] + "ane"
    return f"{prefixes}{'-' if prefixes else ''}{name}"
