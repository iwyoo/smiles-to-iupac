"""Hapto ligands beyond a lone carbocyclic ring (P-69.2.4, P-69.2.6):
acyclic polyenyls and partially bound chains, partially bound or radical
carbocycles, mono-heteroatom five/six-membered rings (pyrrolyl, thiophene),
and fused carbocycles (naphthalene, indenyl, fluorenyl, azulene ...). Names
are built from the ligand's own numbering so that locant sets such as
'(1-3-eta)' can be cited; the metal binds the donors as a delocalised unit.
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._metallacycle import _ring_stem
from ._metallafused import _TEMPLATES, _mappings, _template
from ._numerals import alkane_name
from ._prefix_groups import PrefixNamer
from ._substituents import format_substituent_prefixes

ETA = "η"
DASH = "–"
_VALENCE = {5: 3, 6: 4, 7: 3, 8: 2, 14: 4, 15: 3, 16: 2, 32: 4, 33: 3, 34: 2, 51: 3, 52: 2}
_FUSED = dict(_TEMPLATES)
_FUSED["azulene"] = _template(
    ["1", "2", "3", "3a", "4", "5", "6", "7", "8", "8a"],
    [("1", "2"), ("2", "3"), ("3", "3a"), ("3a", "4"), ("4", "5"), ("5", "6"), ("6", "7"), ("7", "8"),
     ("8", "8a"), ("8a", "1"), ("3a", "8a")],
)
_FUSION_ATOMS = {"3a", "4a", "4b", "7a", "8a", "9a", "10a"}
_HETERO5 = {
    8: "furan", 16: "thiophene", 34: "selenophene", 52: "tellurophene", 7: "1H-pyrrole", 15: "1H-phosphole",
    33: "1H-arsole", 51: "1H-stibole", 5: "1H-borole", 14: "1H-silole", 32: "1H-germole",
}
_HETERO6 = {7: "pyridine", 15: "phosphinine", 33: "arsinine", 51: "stibinine", 5: "borinine"}
_PARENT_WORD = {"indene": "indene", "fluorene": "fluorene", "naphthalene": "naphthalene",
                "anthracene": "anthracene", "phenanthrene": "phenanthrene", "azulene": "azulene"}
_INDICATED = {"indene", "fluorene"}


def lkey(locant):
    m = re.match(r"(\d+)(.*)", str(locant))
    return int(m.group(1)), m.group(2)


def locant_text(locants):
    locs = sorted(locants, key=lkey)
    out, i = [], 0
    while i < len(locs):
        j = i
        while (
            j + 1 < len(locs)
            and str(locs[j]).isdigit()
            and str(locs[j + 1]).isdigit()
            and int(locs[j + 1]) == int(locs[j]) + 1
        ):
            j += 1
        if j - i >= 2:
            out.append(f"{locs[i]}{DASH}{locs[j]}")
        else:
            out.extend(str(x) for x in locs[i : j + 1])
        i = j + 1
    return ",".join(out)


def _kekule(mol):
    k = Chem.Mol(mol)
    Chem.Kekulize(k, clearAromaticFlags=True)
    return k


def valence_gaps(mol, atoms, sigma_orders=None):
    """Free valences per ligand atom; a carbon bonded to the metal by a
    sigma/multiple bond (alkylidene) has that many satisfied by the metal."""
    k = _kekule(mol)
    gaps = {}
    for i in atoms:
        a = k.GetAtomWithIdx(i)
        used = sum(b.GetBondTypeAsDouble() for b in a.GetBonds() if b.GetOtherAtomIdx(i) in atoms) + a.GetTotalNumHs()
        gaps[i] = _VALENCE.get(a.GetAtomicNum(), used) - used
        if sigma_orders and i in sigma_orders and a.GetAtomicNum() == 6:
            gaps[i] -= int(sigma_orders[i])
    if any(g < 0 for g in gaps.values()) or sum(gaps.values()) > 1:
        raise UnsupportedStructure("this electron count of a hapto ligand is not supported here")
    return gaps, k


def _grouped(mol, graph, skeleton_order, allowed, namer):
    grouped: dict = {}
    for position, atom in enumerate(skeleton_order, start=1):
        for n in graph[atom]:
            if n in skeleton_order or n not in allowed:
                continue
            name, compound = namer.name(n, atom)
            entry = grouped.setdefault(name, {"locants": [], "compound": compound})
            entry["locants"].append(position)
    return grouped


def _grouped_labels(mol, graph, mapping, allowed, namer):
    grouped: dict = {}
    for atom, label in mapping.items():
        for n in graph[atom]:
            if n in mapping or n not in allowed:
                continue
            name, compound = namer.name(n, atom)
            entry = grouped.setdefault(name, {"locants": [], "compound": compound})
            entry["locants"].append(label)
    return grouped


def _key(grouped):
    return (
        sorted((lkey(l) for e in grouped.values() for l in e["locants"])),
        [[lkey(l) for l in grouped[k]["locants"]] for k in sorted(grouped)],
    )


def _ring_orders(graph, ring_atoms, start):
    ring_atoms = set(ring_atoms)
    orders = []
    for direction in (0, 1):
        order, prev = [start], None
        while len(order) < len(ring_atoms):
            nbrs = [v for v in graph[order[-1]] if v in ring_atoms and v != prev]
            nxt = nbrs[direction] if len(nbrs) > 1 else nbrs[0]
            prev = order[-1]
            order.append(nxt)
        orders.append(order)
    return orders


def _unsat_stem(size, enes, ynes, ring):
    if ring:
        if ynes:
            raise UnsupportedStructure("a ring triple bond is not supported here")
        return _ring_stem(size, enes, False)
    base = alkane_name(size)
    if not enes and not ynes:
        return base
    stem = base[:-3]
    if enes and ynes:
        e = {1: "en", 2: "dien", 3: "trien", 4: "tetraen"}[len(enes)]
        y = {1: "yne", 2: "diyne", 3: "triyne"}[len(ynes)]
        return f"{stem}{'a' if len(enes) > 1 else ''}-{','.join(map(str, enes))}-{e}-{','.join(map(str, ynes))}-{y}"
    locs = ",".join(map(str, enes or ynes))
    word = {1: "ene", 2: "diene", 3: "triene", 4: "tetraene"} if enes else {1: "yne", 2: "diyne", 3: "triyne"}
    count = len(enes or ynes)
    return f"{stem}{'a' if count > 1 else ''}-{locs}-{word[count]}" if size > 2 or count > 1 else f"{stem}{word[count]}"


def _ending(mol, skeleton):
    """'yl' / 'ido' / 'ylium' from the formal charge on the skeleton atoms
    (IR-10.2.5: anions end in ido, cations in ylium, radicals in yl)."""
    charge = sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in skeleton)
    return "ido" if charge < 0 else "ylium" if charge > 0 else "yl"


def _yl(stem, ending="yl", locant=1):
    return (stem[:-1] if stem.endswith("e") else stem) + f"-{locant}-{ending}"


class Result:
    """A hapto ligand with the locant of every donor atom. For substituent-
    group and ionic names (yl, ido, ylium) the eta term sits between the
    prefixes and the parent stem ('pentamethyl-eta5-cyclopentadienyl'); a
    neutral molecule name starts with it ('eta6-hexamethylbenzene')."""

    def __init__(self, prefix, stem, labels, pi, whole=False, neutral=False, kappa="", force=False):
        self.prefix, self.stem, self.labels, self.pi = prefix, stem, labels, pi
        self.whole, self.neutral, self.kappa, self.force = whole, neutral, kappa, force

    @property
    def name(self):
        return f"{self.prefix}{'-' if self.prefix and self.stem[:1].isdigit() else ''}{self.stem}"

    def needs_locants(self):
        if self.whole:
            return False
        return self.force or _needs_locants(list(self.labels.values()), self.pi)


def render(res, donors=None):
    """Label for the whole ligand (donors=None) or for the subset bound to
    one metal of a bridge; a subset always cites its locants."""
    if donors is not None:
        return f"{locant_text([res.labels[d] for d in donors])}-{ETA}"
    locants = res.needs_locants()
    eta = f"({locant_text(res.labels.values())}-{ETA})-" if locants else f"{ETA}{len(res.labels)}-"
    if res.neutral:
        body = f"{eta}{res.prefix}{res.stem}"
    else:
        body = f"{res.prefix}{'-' if res.prefix else ''}{eta}{res.stem}"
    body += res.kappa
    return f"[{body}]" if locants else body


def _needs_locants(donor_positions, pi_bonds, runs_cover_all=True):
    """Locants are cited when a pi bond is left unbound or the donors form
    several separate runs (P-69.2.4)."""
    donors = set(donor_positions)
    if any(not (a in donors and b in donors) for a, b in pi_bonds):
        return True
    ordered = sorted(donors, key=lkey)
    runs = 1 + sum(1 for x, y in zip(ordered, ordered[1:]) if not (str(x).isdigit() and str(y).isdigit() and int(y) == int(x) + 1))
    return runs > 1


# ---- acyclic -------------------------------------------------------------


def _chain_paths(graph, carbons):
    """All simple carbon paths (tuples), longest first."""
    paths = []

    def walk(path):
        extended = False
        for n in graph[path[-1]]:
            if n in carbons and n not in path:
                walk(path + [n])
                extended = True
        if not extended:
            paths.append(path)

    for start in carbons:
        if sum(1 for n in graph[start] if n in carbons) <= 1:
            walk([start])
    return paths


def chain_label(mol, graph, donor_idx, atoms, sigma=frozenset(), sigma_orders=None):
    if any(mol.GetAtomWithIdx(i).IsInRing() for i in atoms):
        return None
    carbons = {i for i in atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 6}
    if not donor_idx <= carbons:
        raise UnsupportedStructure("only carbon donors are supported in an acyclic hapto ligand")
    gaps, k = valence_gaps(mol, atoms, sigma_orders)
    radical = sum(gaps.values()) == 1
    namer = PrefixNamer(mol, graph, sigma)
    n = len(donor_idx)
    sigma_carbons = {i for i in sigma if mol.GetAtomWithIdx(i).GetAtomicNum() == 6}
    best = None
    for path in _chain_paths(graph, carbons):
        if not donor_idx <= set(path) or not sigma_carbons <= set(path):
            continue
        for order in (path, path[::-1]):
            pos = {a: i + 1 for i, a in enumerate(order)}
            dpos = sorted(pos[d] for d in donor_idx)
            if dpos != list(range(dpos[0], dpos[0] + n)):
                continue
            grouped = _grouped(mol, graph, order, atoms, namer)
            if radical:
                outside = [
                    min(pos[b.GetBeginAtomIdx()], pos[b.GetEndAtomIdx()])
                    for b in k.GetBonds()
                    if b.GetBondTypeAsDouble() == 2.0
                    and b.GetBeginAtomIdx() in pos
                    and b.GetEndAtomIdx() in pos
                    and not (dpos[0] <= pos[b.GetBeginAtomIdx()] <= dpos[-1] and dpos[0] <= pos[b.GetEndAtomIdx()] <= dpos[-1])
                ]
                for free, step in ((dpos[0], 1), (dpos[-1], -1)):
                    inside = [
                        min(free + step * (2 * i - 1), free + step * (2 * i)) for i in range(1, (n - 1) // 2 + 1)
                    ]
                    enes = sorted(inside + outside)
                    key = (-len(order), free, sorted(enes), enes, dpos, _key(grouped))
                    if best is None or key < best[0]:
                        best = (key, order, enes, [], grouped, dpos, pos, free)
            else:
                enes = sorted(
                    min(pos[b.GetBeginAtomIdx()], pos[b.GetEndAtomIdx()])
                    for b in k.GetBonds()
                    if b.GetBondTypeAsDouble() == 2.0 and b.GetBeginAtomIdx() in pos and b.GetEndAtomIdx() in pos
                )
                ynes = sorted(
                    min(pos[b.GetBeginAtomIdx()], pos[b.GetEndAtomIdx()])
                    for b in k.GetBonds()
                    if b.GetBondTypeAsDouble() == 3.0 and b.GetBeginAtomIdx() in pos and b.GetEndAtomIdx() in pos
                )
                key = (-len(order), 0, sorted(enes + ynes), enes, dpos, _key(grouped))
                if best is None or key < best[0]:
                    best = (key, order, enes, ynes, grouped, dpos, pos, 0)
    if best is None:
        raise UnsupportedStructure("this acyclic hapto ligand is not supported here")
    _, order, enes, ynes, grouped, dpos, pos, free = best
    stem = _unsat_stem(len(order), enes, ynes, False)
    if radical:
        if len(order) == 3 and n == 3 and not grouped and free == 1:
            return Result("", "allyl", {d: pos[d] for d in donor_idx}, [], whole=True)
        stem = _yl(stem, _ending(mol, order), free)
    pi = [(e, e + 1) for e in enes + ynes]
    extends = dpos[0] != 1 or dpos[-1] != len(order)
    kappa = ""
    if sigma_carbons:
        words = {1: "yl", 2: "ylidene", 3: "ylidyne"}
        for atom in sorted(sigma_carbons, key=lambda a: pos[a]):
            stem += f"-{pos[atom]}-{words[int(sigma_orders[atom])]}"
        kappa = "-\u03ba" + ("%d" % len(sigma_carbons) if len(sigma_carbons) > 1 else "") + ",".join(
            f"C{pos[a]}" for a in sorted(sigma_carbons, key=lambda a: pos[a])
        )
        if not radical:
            raise UnsupportedStructure("an alkylidene sigma donor on a neutral pi ligand is not supported here")
    return Result(
        format_substituent_prefixes(grouped), stem, {d: pos[d] for d in donor_idx}, pi,
        neutral=not radical, force=(radical and extends) or bool(sigma_carbons), kappa=kappa,
    )


# ---- carbocycle partially bound or radical --------------------------------


def ring_label(mol, graph, ring, donor_idx, atoms, sigma=frozenset()):
    ring_atoms = set(ring)
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring_atoms):
        return hetero_ring_label(mol, graph, ring, donor_idx, atoms, sigma)
    gaps, k = valence_gaps(mol, atoms)
    radical = sum(gaps.values()) == 1
    size, n = len(ring), len(donor_idx)
    namer = PrefixNamer(mol, graph, sigma)
    doubles = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in k.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and b.GetBeginAtomIdx() in ring_atoms and b.GetEndAtomIdx() in ring_atoms
    }
    best = None
    for start in ring_atoms:
        for order in _ring_orders(graph, ring_atoms, start):
            pos = {a: i + 1 for i, a in enumerate(order)}
            dpos = sorted(pos[d] for d in donor_idx)
            if radical:
                if n == size:
                    segment_ok = True
                else:
                    segment_ok = dpos == list(range(1, n + 1))
                if not segment_ok:
                    continue
                enes = [2 * i for i in range(1, (n - 1) // 2 + 1)]
                if n != size:
                    enes += [
                        min(pos[a], pos[b])
                        for a, b in map(tuple, doubles)
                        if min(pos[a], pos[b]) > n and abs(pos[a] - pos[b]) == 1
                    ]
                enes = sorted(enes)
            else:
                enes = []
                for pair in doubles:
                    a, b = tuple(pair)
                    if abs(pos[a] - pos[b]) != 1:
                        enes = None
                        break
                    enes.append(min(pos[a], pos[b]))
                if enes is None:
                    continue
                enes.sort()
            grouped = _grouped(mol, graph, order, atoms, namer)
            key = (enes, _key(grouped), dpos)
            if best is None or key < best[0]:
                best = (key, enes, grouped, dpos, pos)
    if best is None:
        raise UnsupportedStructure("this hapto ring is not supported here")
    _, enes, grouped, dpos, pos = best
    benzene = size == 6 and len(enes) == 3 and not radical
    stem = _ring_stem(size, enes, benzene)
    if radical:
        stem = _yl(stem, _ending(mol, ring_atoms))
    pi = [(e, e + 1) for e in enes]
    return Result(format_substituent_prefixes(grouped), stem, {d: pos[d] for d in donor_idx}, pi, whole=n == size, neutral=not radical)


# ---- hetero monocycle ------------------------------------------------------


def hetero_ring_label(mol, graph, ring, donor_idx, atoms, sigma=frozenset()):
    ring_atoms = set(ring)
    hetero = [i for i in ring_atoms if mol.GetAtomWithIdx(i).GetAtomicNum() != 6]
    size, n = len(ring), len(donor_idx)
    if len(hetero) != 1 or size not in (5, 6):
        raise UnsupportedStructure("only a ring with one heteroatom is supported as a hapto ligand")
    h = hetero[0]
    z = mol.GetAtomWithIdx(h).GetAtomicNum()
    table = _HETERO5 if size == 5 else _HETERO6
    if z not in table:
        raise UnsupportedStructure("this ring heteroatom is not supported in a hapto ligand")
    gaps, _ = valence_gaps(mol, atoms)
    radical = sum(gaps.values()) == 1
    if radical and (gaps[h] != 1 or size != 5 or z in (8, 16, 34, 52)):
        raise UnsupportedStructure("this radical hapto ring is not supported here")
    namer = PrefixNamer(mol, graph, sigma)
    best = None
    for order in _ring_orders(graph, ring_atoms, h):
        pos = {a: i + 1 for i, a in enumerate(order)}
        grouped = _grouped(mol, graph, order, atoms, namer)
        key = _key(grouped)
        if best is None or key < best[0]:
            best = (key, grouped, sorted(pos[d] for d in donor_idx), pos)
    _, grouped, dpos, pos = best
    parent = table[z]
    stem = _yl(parent, _ending(mol, ring_atoms)) if radical else parent
    return Result(format_substituent_prefixes(grouped), stem, {d: pos[d] for d in donor_idx}, [], whole=n == size, neutral=not radical)


# ---- fused carbocycles -----------------------------------------------------


def fused_label(mol, graph, system, donor_idx, atoms, sigma=frozenset()):
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in system):
        raise UnsupportedStructure("only a fused carbocycle is supported as a hapto ligand")
    gaps, k = valence_gaps(mol, atoms)
    radical = sum(gaps.values()) == 1
    sub = {a: {x for x in graph[a] if x in system} for a in system}
    namer = PrefixNamer(mol, graph, sigma)
    sp3 = [
        i
        for i in system
        if not any(
            b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtomIdx(i) in system for b in k.GetAtomWithIdx(i).GetBonds()
        )
    ]
    best = None
    for tname, (labels, adj) in _FUSED.items():
        if len(labels) != len(system):
            continue
        for mapping in _mappings(system, sub, labels, adj):
            grouped = _grouped_labels(mol, graph, mapping, atoms, namer)
            dlabels = sorted((mapping[d] for d in donor_idx), key=lkey)
            sp3_labels = sorted((mapping[i] for i in sp3), key=lkey)
            if tname in _INDICATED:
                if radical:
                    cands = [mapping[d] for d in donor_idx if mapping[d] not in _FUSION_ATOMS]
                    if not cands:
                        continue
                    attach = min(cands, key=lkey)
                    indicated = attach
                else:
                    if len(sp3_labels) != 1:
                        continue
                    attach, indicated = None, sp3_labels[0]
            else:
                if radical or sp3_labels:
                    continue
                attach, indicated = None, None
            key = (lkey(indicated) if indicated else (0, ""), _key(grouped), [lkey(l) for l in dlabels])
            if best is None or key < best[0]:
                best = (key, tname, grouped, dlabels, attach, indicated, mapping)
        if best is not None:
            break
    if best is None:
        raise UnsupportedStructure("this fused hapto ligand is not supported here")
    _, tname, grouped, dlabels, attach, indicated, mapping = best
    parent = _PARENT_WORD[tname]
    hydrogen = f"{indicated}H-" if indicated else ""
    if radical:
        stem = f"{hydrogen}{_yl(parent, _ending(mol, system), attach)}"
    else:
        stem = f"{hydrogen}{parent}"
    pi = [
        (mapping[b.GetBeginAtomIdx()], mapping[b.GetEndAtomIdx()])
        for b in k.GetBonds()
        if b.GetBondTypeAsDouble() >= 2.0 and b.GetBeginAtomIdx() in mapping and b.GetEndAtomIdx() in mapping
    ]
    return Result(format_substituent_prefixes(grouped), stem, {d: mapping[d] for d in donor_idx}, pi, neutral=not radical)
