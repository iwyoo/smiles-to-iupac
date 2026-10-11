"""Organic rings and ring systems ortho fused to a fullerene or modified fullerene, cited as fusion prefixes with primed
locants in front of the fullerene parent (P-27.6.1). The fusion bond stays single and the fusion atoms take no exo double
bond, so the ring component is judged without that bond: its remaining saturated positions take indicated hydrogen
when odd in number and hydro prefixes otherwise (P-31.1.4.2.4)."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fused_numbering import _locant_key
from ._fusion_components import _RETAINED_PREFIX, identify, skeleton
from ._numerals import multiplying_prefix

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


@dataclass
class Placed:
    ring: Ring
    cage: tuple
    locants: tuple
    option: tuple
    mark: str


def ring_component(mol, kekule, graph, a, b, atoms):
    """The component fused at the bond a-b, with the numberings that give its fusion atoms and its indicated hydrogen
    the lowest locants; None when it cannot be named as one fusion component."""
    members = sorted(atoms | {a, b})
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
    try:
        components = identify(sub)
    except UnsupportedStructure:
        return None
    if not components or len({c.bracketed_prefix() for c in components}) > 1:
        return None
    component = components[0]

    def saturated(x):
        atom = mol.GetAtomWithIdx(x)
        return (
            atom.GetSymbol() not in _DONORS
            and atom.GetTotalNumHs() > 0
            and all(kekule.GetBondBetweenAtoms(x, n).GetBondType() != Chem.BondType.DOUBLE for n in graph[x] if n in index)
        )

    options = []
    for numbering in component.numberings:
        locant = {x: numbering[index[x]] for x in members}
        positions = sorted((locant[x] for x in atoms if saturated(x)), key=_locant_key)
        odd = len(positions) % 2
        indicated, hydro = tuple(positions[:odd]), tuple(positions[odd:])
        key = (
            tuple(sorted(_locant_key(locant[x]) for x in (a, b))),
            tuple(_locant_key(p) for p in indicated),
            tuple(_locant_key(p) for p in hydro),
        )
        options.append((key, (locant[a], locant[b], indicated, hydro, locant)))
    best = min(key for key, _ in options)
    kind = component.kind in ("mono", "attached_only")
    return Ring(
        a=a,
        b=b,
        atoms=frozenset(atoms),
        prefix=component.prefix,
        bracketed=component.bracketed_prefix(),
        cite=not (kind and all(e == "C" for e in component.elements)),
        compound=bool(component.hetero_text)
        or (kind and component.name not in _RETAINED_PREFIX and any(e != "C" for e in component.elements)),
        options=[option for key, option in options if key == best],
    )


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
    return hydro_text, indicated, "".join(texts)
