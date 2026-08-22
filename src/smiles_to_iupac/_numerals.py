"""Numerical terms for chain lengths (IUPAC 2013 Recommendations, "the Blue Book",
P-14.2.1, Table 1.4; https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf).
"""

_UNITS = ["", "mono", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona"]
_UNITS_COMPOSITE = ["", "hen", "do", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona"]
_TENS = ["", "deca", "icosa", "triaconta", "tetraconta", "pentaconta",
         "hexaconta", "heptaconta", "octaconta", "nonaconta"]
_HUNDREDS = ["", "hecta", "dicta", "tricta", "tetracta", "pentacta",
             "hexacta", "heptacta", "octacta", "nonacta"]
_THOUSANDS = ["", "kilia", "dilia", "trilia", "tetralia", "pentalia",
              "hexalia", "heptalia", "octalia", "nonalia"]

_RETAINED_ALKANES = {1: "methane", 2: "ethane", 3: "propane", 4: "butane"}
_RETAINED_ALKYLS = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}


def numerical_term(n: int) -> str:
    """Basic numerical term (multiplying prefix) for n, per P-14.2.1 / Table 1.4."""
    if n <= 0 or n > 9999:
        raise ValueError(f"numerical term not defined for n={n}")
    if n < 10:
        return _UNITS[n]
    if n == 11:
        return "undeca"

    units, tens, hundreds, thousands = n % 10, (n // 10) % 10, (n // 100) % 10, (n // 1000) % 10
    composite = tens or hundreds or thousands

    parts = []
    if units:
        parts.append(_UNITS_COMPOSITE[units] if composite else _UNITS[units])
    if tens == 1 and units == 1:
        # P-14.2.1.1.1: 11 is the irregular 'undeca', not 'hen' + 'deca'.
        parts = ["undeca"]
    elif tens:
        term = _TENS[tens]
        # P-14.2.1.2: the leading 'i' of 'icosa' elides after a vowel.
        if term == "icosa" and parts and parts[-1][-1] in "aeiou":
            term = term[1:]
        parts.append(term)
    if hundreds:
        parts.append(_HUNDREDS[hundreds])
    if thousands:
        parts.append(_THOUSANDS[thousands])
    return "".join(parts)


def alkane_name(n: int) -> str:
    """Name of the unbranched acyclic hydrocarbon with n carbons (P-21.2.1)."""
    if n in _RETAINED_ALKANES:
        return _RETAINED_ALKANES[n]
    term = numerical_term(n)
    # P-21.2.1: numerical term + 'ane', eliding the term's terminal 'a'.
    return term[:-1] + "ane" if term.endswith("a") else term + "ane"


def alkyl_name(n: int) -> str:
    """Name of the substituent group -[CH2]n-1-CH3 (P-29.3.2.1)."""
    if n in _RETAINED_ALKYLS:
        return _RETAINED_ALKYLS[n]
    term = numerical_term(n)
    return term[:-1] + "yl" if term.endswith("a") else term + "yl"


_KIS_IRREGULAR = {2: "bis", 3: "tris"}


def multiplying_prefix(n: int, compound: bool = False) -> str:
    """Multiplying prefix for n identical substituent prefixes: the basic
    numerical term (P-14.2.1) for simple substituents, or the irregular
    'bis'/'tris' or regular '...kis' series (P-14.2.2) for compound
    substituents, avoiding ambiguity with a substituent's own internal
    multiplying prefixes."""
    if not compound:
        return numerical_term(n)
    if n in _KIS_IRREGULAR:
        return _KIS_IRREGULAR[n]
    return numerical_term(n) + "kis"
