"""Unsaturation, parent text and skeletal-modification prefixes of natural-product names (P-101.3, P-101.6)."""

from ._common import multiplied_word
from ._np_core import loc_key
from ._numerals import multiplying_prefix


_RETAINED_ENDINGS = {"morphinan", "ibogamine"}


def stem_info(name):
    if name in _RETAINED_ENDINGS:
        return name, None
    for ending, kind in (("anine", "anine"), ("ane", "ane"), ("an", "an")):
        if name.endswith(ending):
            return name[: -len(ending)], kind
    return name, None


def unsaturation(cand, view):
    """([ene bonds (a, b, order)], [hydro locants], [dehydro locants]) relative to the (modified) parent (P-101.6)."""
    skel, mapping = cand.skel, cand.mapping
    parent = skel.parent
    aromatic_m = {i for bond in view.aromatic_bonds for i in bond}
    aromatic = {loc for loc, atom in mapping.items() if loc in parent.aromatic_atoms or atom in aromatic_m}
    _, kind = stem_info(parent.name)
    saturated = {
        loc: all(skel.order.get(frozenset((loc, n)), 1) == 1 for n in skel.adj[loc]) and loc not in parent.aromatic_atoms
        for loc in mapping
    }
    ene, hydro, dehydro = [], [], []
    covered = set()
    if kind:
        ring = {loc for loc in aromatic if loc not in parent.aromatic_atoms and saturated[loc]}
        pairs = [
            tuple(bond) for bond in view_aromatic_pairs(cand, view)
            if all(x in ring for x in bond)
        ]
        chosen = _kekule(ring, pairs, sorted(skel.adj, key=lambda x: loc_key(x)))
        if chosen is not None:
            ene += [(a, b, 2) for a, b in chosen]
            covered = set(ring)
    for bond, order in skel.order.items():
        a, b = tuple(bond)
        if a in aromatic and b in aromatic:
            continue
        other = view.order.get(frozenset((mapping[a], mapping[b])))
        if other is None or other == order:
            continue
        if other > order:
            if kind and order == 1 and saturated[a] and saturated[b]:
                ene.append((a, b, other))
            else:
                dehydro += [a, b] * (other - order)
        else:
            hydro += [a, b] * (order - other)
    for loc in aromatic:
        if loc in covered or (loc in parent.aromatic_atoms and mapping[loc] in aromatic_m):
            continue
        d_p = 1 if loc in parent.aromatic_atoms else sum(skel.order.get(frozenset((loc, n)), 1) - 1 for n in skel.adj[loc])
        d_m = sum(
            view.order[frozenset((mapping[loc], mapping[n]))] - 1
            for n in skel.adj[loc]
            if frozenset((mapping[loc], mapping[n])) in view.order
        )
        if d_m < d_p:
            hydro += [loc] * (d_p - d_m)
        elif d_m > d_p:
            dehydro += [loc] * (d_m - d_p)
    return ene, sorted(hydro, key=loc_key), sorted(dehydro, key=loc_key)


def locant_pair(skel, final, a, b):
    ordered = sorted(skel.adj, key=lambda x: loc_key(final(x)))
    low, high = sorted((a, b), key=lambda x: loc_key(final(x)))
    if abs(ordered.index(low) - ordered.index(high)) == 1:
        return final(low)
    return f"{final(low)}({final(high)})"


def core_text(parent_name, enes, ynes, suffix, locants):
    """Parent name with ene/yne endings and the principal suffix (P-31.1.4, P-101.6.1)."""
    stem, kind = stem_info(parent_name)
    if not kind:
        text = stem
        if suffix:
            if suffix[0] in "aeiouy" and text.endswith("e"):
                text = text[:-1]
            text += f"-{','.join(locants)}-{suffix}"
        return text
    segments = []
    if enes:
        segments.append((enes, multiplied_word(len(enes), {"anine": "enine"}.get(kind, "ene"))))
    if ynes:
        segments.append((ynes, multiplied_word(len(ynes), {"anine": "ynine"}.get(kind, "yne"))))
    if not segments:
        segments.append(([], kind))
    if suffix:
        segments.append((locants, suffix))
    first_words = segments[0][1]
    text = stem + ("a" if segments[0][0] and first_words.startswith(("di", "tri", "tetra", "penta", "hexa")) else "")
    for i, (locs, word) in enumerate(segments):
        if word.endswith("e") and i + 1 < len(segments) and segments[i + 1][1][0] in "aeiouy":
            word = word[:-1]
        text += (f"-{','.join(locs)}-" if locs else "") + word
    return text


def _multiplied(items, word):
    if not items:
        return None
    mult = "" if len(items) == 1 else multiplying_prefix(len(items), compound=False)
    return f"{mult}{word}"


def op_prefixes(cand, final, cyclo=None):
    """The nondetachable prefixes of the skeletal modifications in citation order (P-101.3.7.2)."""
    ops = cand.skel.ops
    parts = []
    cyclo = cyclo or []
    if cyclo:
        parts.append(f"{':'.join(cyclo)}-{_multiplied(cyclo, 'cyclo')}")
    secos = [f"{final(op[1])},{final(op[2])}" for op in ops if op[0] == "seco"]
    if secos:
        parts.append(f"{':'.join(secos)}-{_multiplied(secos, 'seco')}")
    homo = sorted((final(op[3]) for op in ops if op[0] == "homo"), key=loc_key)
    if homo:
        parts.append(f"{','.join(homo)}-{_multiplied(homo, 'homo')}")
    nor = sorted((final(op[1]) for op in ops if op[0] == "nor"), key=loc_key)
    if nor:
        parts.append(f"{','.join(nor)}-{_multiplied(nor, 'nor')}")
    return parts


def view_aromatic_pairs(cand, view):
    image = {a: loc for loc, a in cand.mapping.items()}
    return [
        frozenset(image[i] for i in bond)
        for bond in view.aromatic_bonds
        if all(i in image for i in bond)
    ]


def _kekule(atoms, pairs, ordering):
    """Double bonds of the Kekulé structure with the lowest locants, simple before compound (P-31.1.4.2.4)."""
    if not atoms:
        return None
    order = sorted(atoms, key=loc_key)
    best = []

    def solve(remaining, chosen):
        if not remaining:
            best.append(list(chosen))
            return
        first = remaining[0]
        for a, b in pairs:
            other = b if a == first else a if b == first else None
            if other is not None and other in remaining:
                rest = [x for x in remaining if x not in (first, other)]
                chosen.append((first, other))
                solve(rest, chosen)
                chosen.pop()

    solve(order, [])
    if not best:
        return None
    position = {a: i for i, a in enumerate(ordering)}

    def citation(bond):
        low, high = sorted(bond, key=loc_key)
        if abs(position[low] - position[high]) == 1:
            return (loc_key(low),)
        return (loc_key(low), loc_key(high))

    return min(best, key=lambda bonds: sorted(citation(b) for b in bonds))
