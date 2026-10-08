"""Retained formic, acetic and oxalic acid names and the anion and acyl names built on them (P-65.1.1.1, P-65.1.7.2.1,
P-65.1.8, P-15.1.8.2.1). Acetic acid keeps its name under substitution, formic acid only for groups P-65.1.8.1 leaves
to carbonic acid names, oxalic acid never; the single position needs no locant and later prefixes are enclosed.
"""

import re

from ._common import alpha_sort_key
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._substituents import prefix_multiplier

_SUBSTITUTED_RETAINED_ACYL = re.compile(r"(?<=.)(?:acetyl|formyl)$")
_HYDRIDE_ACYL = re.compile(r"(?<=.)(?:thioyl|selenoyl|telluroyl|imidoyl|hydrazonoyl|carbonyl|sulfonyl|sulfinyl|selenonyl|seleninyl|telluronyl|tellurinyl)$")
_ENDINGS = {"acid": "ic acid", "anion": "ate", "acyl": "yl", "acylium": "ylium"}
_OXALIC = {"acid": "oxalic acid", "anion": "oxalate", "acyl": "oxalyl", "acylium": "oxalylium"}
_OXAMIC = {"acid": "oxamic acid", "anion": "oxamate", "acyl": "oxamoyl", "acylium": "oxamoylium"}
# P-65.1.8.1: these replace the hydrogen of formic acid only in carbonic acid names; hydroxy, alkoxy and
# alkylsulfanyl groups make a carbonic acid ester in the same way
_CARBONIC_PREFIXES = {
    "hydroperoxy", "sulfanyl", "selanyl", "tellanyl", "fluoro", "chloro", "bromo", "iodo", "azido", "isocyano",
    "cyano", "isocyanato", "isothiocyanato", "isoselenocyanato", "isotellurocyanato", "amino", "hydrazinyl", "hydroxy",
}
_CARBONIC_ENDINGS = ("oxy", "sulfanyl", "selanyl", "tellanyl", "amino", "anilino", "peroxy", "azanyl")


def _formic_substitutable(names):
    return not any(n in _CARBONIC_PREFIXES or n.endswith(_CARBONIC_ENDINGS) for n in names)


def single_site_prefixes(grouped):
    """Prefixes of a parent with one substitutable position (P-16.5.1.3.1): no
    locants; the first prefix is enclosed only when it is a compound or
    locant-bearing name, every later one always, multipliers outside."""
    parts = []
    for position, name in enumerate(sorted(grouped, key=alpha_sort_key)):
        info = grouped[name]
        count = len(info["locants"])
        compound = info["compound"]
        multiplier, enclosed = prefix_multiplier(count, name, compound) if count > 1 else ("", False)
        parts.append(multiplier + (enclose(name) if compound or position or enclosed else name))
    return "".join(parts)


def retained_chain_acid(grouped, chain_length, ene_locants, yne_locants, count, family):
    """Name of an acid, anion or acyl group on a one- or two-carbon chain with
    a retained name, or None; `family` is 'acid', 'anion', 'acyl' or 'acylium'."""
    if ene_locants or yne_locants:
        return None
    if chain_length == 2 and count == 2 and not grouped:
        return _OXALIC[family]
    if (
        chain_length == 2
        and count == 1
        and set(grouped) == {"amino", "oxo"}
        and all(info["locants"] == [2] and not info["compound"] for info in grouped.values())
    ):
        return _OXAMIC[family]
    if chain_length == 2 and count == 1:
        stem = "acet"
    elif chain_length == 1 and count == 1 and _formic_substitutable(grouped):
        stem = "form"
    else:
        return None
    return single_site_prefixes(grouped) + stem + _ENDINGS[family]


def is_compound_acyl(name):
    """Whether an acyl prefix name needs enclosing marks: a substituted group, or one that contains the
    name of a parent hydride (P-16.5.1.1, P-16.5.1.4)."""
    return (
        any(ch.isdigit() for ch in name)
        or "(" in name
        or "[" in name
        or bool(_SUBSTITUTED_RETAINED_ACYL.search(name))
        or bool(_HYDRIDE_ACYL.search(name))
    )
