"""Organic rings and ring systems ortho fused to a fullerene or modified fullerene, cited as fusion prefixes with primed
locants in front of the fullerene parent (P-27.6.1). The fusion bond stays single and the fusion atoms take no exo double
bond, so the ring component is judged without that bond: its remaining saturated positions take indicated hydrogen
when odd in number and hydro prefixes otherwise (P-31.1.4.2.4)."""

import re
from dataclasses import dataclass

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents, unsaturation_suffix
from ._fused_numbering import HETERO_RANK, _locant_key
from ._fusion_components import _RETAINED_PREFIX, identify, skeleton
from ._fusion_name import fused_prefix_data
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._substituents import name_branch, prefix_multiplier

_PRIMES = ("", "′", "″", "‴", "⁗")
_DONORS = frozenset({"O", "S", "Se", "Te"})


@dataclass
class Ring:
    """A fused component: `options` are (locant of a, locant of b, indicated hydrogen, hydro, {atom: locant}) of the
    numberings that tie."""

    a: int
    b: int
    atoms: frozenset
    prefix: str
    bracketed: str
    cite: bool
    compound: bool
    options: list
    branches: tuple = ()
    derived: bool = False


@dataclass
class Placed:
    ring: Ring
    cage: tuple
    locants: tuple
    option: tuple
    mark: str


def ring_component(mol, kekule, graph, a, b, atoms):
    """The component fused at the bond a-b, with the numberings that give its fusion atoms, its indicated hydrogen and its
    substituents the lowest locants; None when it cannot be named as one fusion component. The acyclic atoms of the
    group are substituents of the component, cited with its locants."""
    ring_atoms = _fused_ring_atoms(mol, atoms, a, b)
    members = sorted(ring_atoms | {a, b})
    index = {x: i for i, x in enumerate(members)}
    sub = skeleton(
        [mol.GetAtomWithIdx(x).GetSymbol() for x in members],
        [(index[x], index[y]) for x, y in graph.edges if x in index and y in index],
    )
    ring_info = sub.GetRingInfo()
    if any(not ring_info.NumAtomRings(i) for i in range(len(members))):
        return None
    if any(kekule.GetBondBetweenAtoms(x, n).GetBondType() == Chem.BondType.DOUBLE for x in (a, b) for n in graph[x] if n in atoms):
        return None
    branches = _branches(mol, kekule, graph, ring_atoms, atoms - ring_atoms)
    if branches is None:
        return None
    try:
        components = identify(sub)
    except UnsupportedStructure:
        components = []
    leading, texts = None, None
    if components:
        if len({c.bracketed_prefix() for c in components}) > 1:
            return None
        component = components[0]
        numberings, prefix, bracketed = component.numberings, component.prefix, component.bracketed_prefix()
        cite = not (component.kind in ("mono", "attached_only") and all(e == "C" for e in component.elements))
        compound = bool(component.hetero_text) or (
            component.kind in ("mono", "attached_only")
            and component.name not in _RETAINED_PREFIX
            and any(e != "C" for e in component.elements)
        )
    else:
        try:
            bracketed, numberings = fused_prefix_data(sub)
        except UnsupportedStructure:
            bridged = _bridged_alicyclic(mol, kekule, sub, members, index)
            if bridged is None:
                return None
            numberings, leading, texts = bridged
            bracketed = texts[0]
        prefix, cite, compound = re.sub(r"\[[^\]]*\]", "", bracketed), True, True
    carrying = {x for x, _, _ in branches}
    if any(mol.GetAtomWithIdx(x).GetSymbol() in _DONORS for x in carrying):
        return None

    def saturated(x):
        atom = mol.GetAtomWithIdx(x)
        return (
            atom.GetSymbol() not in _DONORS
            and (atom.GetTotalNumHs() > 0 or x in carrying)
            and all(kekule.GetBondBetweenAtoms(x, n).GetBondType() != Chem.BondType.DOUBLE for n in graph[x] if n in index)
        )

    options = []
    for number, numbering in enumerate(numberings):
        locant = {x: numbering[index[x]] for x in members}
        positions = [] if texts else sorted((locant[x] for x in ring_atoms if saturated(x)), key=_locant_key)
        odd = len(positions) % 2
        indicated, hydro = tuple(positions[:odd]), tuple(positions[odd:])
        key = (
            leading[number] if texts else (),
            tuple(sorted(_locant_key(locant[x]) for x in (a, b))),
            tuple(_locant_key(p) for p in indicated),
            tuple(_locant_key(p) for p in hydro),
            tuple(sorted(_locant_key(locant[x]) for x, _, _ in branches)),
            tuple(sorted((alpha_sort_key(name), _locant_key(locant[x])) for x, name, _ in branches)),
        )
        options.append((key, (locant[a], locant[b], indicated, hydro, locant)))
    best = min(key for key, _ in options)
    if texts:
        bracketed = texts[next(number for number, (key, _) in enumerate(options) if key == best)]
        prefix = re.sub(r"[^A-Za-z]", "", bracketed)
    return Ring(
        a=a,
        b=b,
        atoms=frozenset(atoms),
        prefix=prefix,
        bracketed=bracketed,
        cite=cite,
        compound=compound,
        options=[option for key, option in options if key == best],
        branches=tuple(branches),
        derived=bool(texts),
    )


def _bridged_alicyclic(mol, kekule, sub, members, index):
    """(numberings, ranking keys, prefix texts) of a bridged alicyclic ring system: its von Baeyer name with 'ane' changed
    to 'ano', heteroatoms as 'a' prefixes, in its own numbering with the heteroatoms lowest (P-27.6.1, P-23, P-25.3.1.3);
    None when it has none."""
    from ._spiro_union import _replacement_prefixes, _von_baeyer, _von_baeyer_ene

    if sub.GetRingInfo().NumRings() < 2:
        return None
    try:
        stem, numberings = _von_baeyer(sub)
        hetero = [x for x in members if mol.GetAtomWithIdx(x).GetAtomicNum() != 6]
        keys, texts = [], []
        for numbering in numberings:
            locant = {x: numbering[index[x]] for x in members}
            ene = _von_baeyer_ene(kekule, sub, members, numbering)
            if ene:
                body, needs_a = unsaturation_suffix(ene, ())
                name = f"{stem[:-3]}{'a' if needs_a else ''}-{body}"
            else:
                name = stem
            replaced = _replacement_prefixes(mol, hetero, locant) if hetero else ""
            ranked = sorted((HETERO_RANK[mol.GetAtomWithIdx(x).GetSymbol()], _locant_key(locant[x])) for x in hetero)
            keys.append((tuple(sorted(k for _, k in ranked)), tuple(k for _, k in ranked), tuple(ene)))
            texts.append(replaced + name[:-1] + "o")
    except UnsupportedStructure:
        return None
    return numberings, keys, texts


def _fused_ring_atoms(mol, atoms, a, b):
    """The atoms of the ring on the fusion bond and of the rings fused to it; a ring joined to them by a single bond or
    a spiro atom belongs to a substituent."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= atoms | {a, b}]
    system = [r for r in rings if a in r and b in r]
    grown = True
    while grown:
        grown = False
        for r in rings:
            if r not in system and any(len(r & s) >= 2 for s in system):
                system.append(r)
                grown = True
    return set().union(*system) - {a, b} if system else set()


def _branches(mol, kekule, graph, ring_atoms, loose):
    """(ring atom, name, is compound) of every acyclic group on the ring atoms, None when a group cannot be cited."""
    if not loose:
        return []
    adj = adjacency(mol)
    halogens = halogen_substituents(mol)
    found, covered = [], set()
    for x in sorted(ring_atoms):
        for n in graph[x]:
            if n in loose:
                if kekule.GetBondBetweenAtoms(x, n).GetBondType() != Chem.BondType.SINGLE:
                    return None
                try:
                    name, compound = name_branch(adj, n, x, halogens, frozenset(), mol=mol)
                except UnsupportedStructure:
                    return None
                found.append((x, name, compound))
                stack = [n]
                while stack:
                    node = stack.pop()
                    if node not in covered:
                        covered.add(node)
                        stack.extend(m for m in graph[node] if m != x and m in loose)
    if covered != set(loose) or any(
        mol.GetAtomWithIdx(i).GetIsotope() or mol.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for i in loose
    ):
        return None
    return found


def layout(rings, numbering, edit):
    """The fused rings of an embedding grouped by prefix in alphanumerical order, each group by its lower fullerene
    locant, primed in that order; every ring with the locants of its fusion atoms cited in the direction of the
    fullerene locants."""
    groups = {}
    for atom_a, atom_b, (kind, _, i) in edit.fusions:
        if kind:
            continue
        ring = rings[i]
        ascending = numbering[atom_a] < numbering[atom_b]
        cage = tuple(sorted((numbering[atom_a], numbering[atom_b])))

        def cited(option):
            first, second = (option[0], option[1]) if ascending else (option[1], option[0])
            return _locant_key(first), _locant_key(second)

        option = min(ring.options, key=cited)
        locants = (option[0], option[1]) if ascending else (option[1], option[0])
        groups.setdefault((ring.prefix, ring.bracketed), []).append([ring, cage, locants, option])
    ordered = [sorted(groups[key], key=lambda item: item[1]) for key in sorted(groups)]
    prime = 0
    for group in ordered:
        for item in group:
            prime += 1
            item.append(_PRIMES[prime] if prime < len(_PRIMES) else "′" * prime)
    return [[Placed(*item) for item in group] for group in ordered]


def ring_key(groups):
    return tuple(_locant_key(l) for group in groups for item in group for l in item.locants)


def ring_positions(groups):
    """The primed locant of every atom of the fused rings."""
    return {atom: f"{locant}{item.mark}" for group in groups for item in group for atom, locant in item.option[4].items() if atom not in (item.ring.a, item.ring.b)}


def substituent_text(groups):
    """The substituent prefixes of the placed rings with their primed locants, in alphanumerical order (P-14.5)."""
    cited = {}
    for members in groups:
        for item in members:
            for atom, name, compound in item.ring.branches:
                entry = cited.setdefault(name, [compound, []])
                entry[1].append(f"{item.option[4][atom]}{item.mark}")
    parts = []
    for name in sorted(cited, key=alpha_sort_key):
        compound, locants = cited[name]
        locants = sorted(locants, key=lambda l: _locant_key(l.rstrip("′″‴⁗")))
        if len(locants) == 1:
            text = enclose(name) if compound else name
        else:
            multiplier, enclosed = prefix_multiplier(len(locants), name, compound)
            text = multiplier + (enclose(name) if enclosed else name)
        parts.append(f"{','.join(locants)}-{text}")
    return "-".join(parts)


def ring_text(groups):
    """(hydro prefixes, indicated hydrogen entries, fusion prefixes) of the placed rings."""
    texts, indicated, hydro = [], [], []
    for members in groups:
        sets = []
        for item in members:
            if item.ring.cite:
                sets.append(f"{item.locants[0]}{item.mark},{item.locants[1]}{item.mark}:{item.cage[0]},{item.cage[1]}")
            else:
                sets.append(f"{item.cage[0]},{item.cage[1]}")
            prime = _PRIMES.index(item.mark) if item.mark in _PRIMES else len(item.mark)
            indicated += [((_locant_key(p)[0], 0, prime), f"{p}{item.mark}H") for p in item.option[2]]
            hydro += [((_locant_key(p), prime), f"{p}{item.mark}") for p in item.option[3]]
        ring = members[0].ring
        name = ring.bracketed
        if len(members) > 1:
            name = multiplying_prefix(len(members), compound=True) + f"({name})" if ring.compound else multiplying_prefix(len(members)) + name
        texts.append(f"{name}[{(';' if ring.cite else ':').join(sets)}]")
    hydro_text = ""
    if hydro:
        locants = [text for _, text in sorted(hydro)]
        hydro_text = ",".join(locants) + "-" + (multiplying_prefix(len(locants)) if len(locants) > 1 else "") + "hydro"
    return hydro_text, indicated, "".join(("-" if i and text[0].isdigit() else "") + text for i, text in enumerate(texts))
