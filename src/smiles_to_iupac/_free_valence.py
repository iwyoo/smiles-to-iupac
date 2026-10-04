"""Free-valence suffixes of substituent groups (P-29.2): 'yl', 'ylidene' and 'ylidyne' mark single, double and
triple attachment; they are cited in that order, each multiplied for several valences of one kind."""

from ._common import multiplied_word

SUFFIX_OF_ORDER = {1: "yl", 2: "ylidene", 3: "ylidyne"}
_VOWELS = "aeiouy"


def suffix_of(order):
    """'yl' / 'ylidene' / 'ylidyne' for a bond order (1, 2, 3 or the float equivalents); None for any other order."""
    return SUFFIX_OF_ORDER.get(int(order)) if order == int(order) else None


def valence_word(count, order=1):
    """'yl', 'diyl', 'triyl', 'ylidene', 'diylidene', ... for `count` valences of one bond order."""
    return multiplied_word(count, SUFFIX_OF_ORDER[order])


def citation(by_order):
    """[(locants text, word)] in the order yl, ylidene, ylidyne from {order: [locant, ...]}."""
    return [(",".join(str(loc) for loc in locs), valence_word(len(locs), order)) for order, locs in sorted(by_order.items())]


def attach(stem, by_order, elide=True):
    """`stem` with its free valences cited: 'cyclohexan-1-yl-2-ylidene', 'propane-1,3-diyl'."""
    pieces = citation(by_order)
    if elide and stem.endswith("e") and pieces[0][1][0] in _VOWELS:
        stem = stem[:-1]
    return stem + "".join(f"-{locs}-{word}" for locs, word in pieces)
