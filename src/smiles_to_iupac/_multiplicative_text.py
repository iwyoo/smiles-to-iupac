"""Text helpers for multiplicative names (P-16.5 enclosing marks, P-16.3
multiplying prefixes) shared by the unit, linker and assembly modules.
"""

import re

from ._common import superscript_locant
from ._numerals import numerical_term

_OPENERS = {"(": 0, "[": 1, "{": 2}
_MARKS = ["()", "[]", "{}"]


_ASSEMBLY_BRACKETS = re.compile(
    r"\[\d+,\d+'-bi(?:\([a-z]+\)|[a-z]+)\]"  # [1,1'-biphenyl]
    r"|\[[\d,:⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ]+-(?:ter|quater|quinque|sexi)[a-z]+\]"  # [1¹,2¹:2⁴,3¹-terphenyl]
    r"|\[[\d.^,':\u2032]+\]"  # von Baeyer, spiro and primed fusion locants: [3.3.1.1^3,7], [4.5], [3',4':6,7]
    r"|\[[\d,']*-?[a-z]{1,2}(?:,[a-z]{1,2})*\]"  # fusion descriptors: [b], [3,2-b], [b,f]
)


_ADDED_HYDROGEN = re.compile(r"\(\d+[a-z]?H\)")


def nesting_text(text):
    """`text` without the marks that P-16.5.4.1 leaves out of the nesting order: the descriptor brackets of ring
    assembly, von Baeyer, spiro and fusion names, and the parentheses of added hydrogen."""
    return _ADDED_HYDROGEN.sub("", _ASSEMBLY_BRACKETS.sub("", text))


def enclose(text):
    depth = deepest = 0
    for ch in nesting_text(text):
        if ch in _OPENERS:
            depth += 1
            deepest = max(deepest, depth)
        elif ch in ")]}":
            depth -= 1
    level = deepest % 3
    left, right = _MARKS[level]
    return f"{left}{text}{right}"


def multiplier_word(count, use_bis):
    if count == 1:
        return ""
    if use_bis:
        return {2: "bis", 3: "tris"}.get(count, numerical_term(count) + "kis")
    return numerical_term(count)


def primed_locants(locant, count):
    return ",".join(str(locant) + "'" * i for i in range(count))


_NUMERIC_STEMS = ("dec", "undec", "dodec", "tridec", "tetradec", "pentadec", "hexadec", "heptadec", "octadec", "nonadec", "icos")


def needs_enclosing(text):
    return any(ch.isdigit() or ch == "-" for ch in text) or text.startswith(_NUMERIC_STEMS)


def unit_phrase(text, tail):
    if needs_enclosing(text):
        return f"({text}{tail})" if tail.startswith(" ") else f"({text}){tail}"
    return text + tail


class PrimedLocant:
    """A locant of a ring assembly: number plus prime count (2 -> 2'), ordered 1 < 1' < 2."""

    def __init__(self, primes, number):
        self.primes, self.number = primes, number

    def _key(self):
        return float(self.number), self.primes

    def __eq__(self, other):
        return isinstance(other, PrimedLocant) and self._key() == other._key()

    def __lt__(self, other):
        return self._key() < other._key()

    def __hash__(self):
        return hash(self._key())

    def __str__(self):
        return f"{self.number}{chr(39) * self.primes}"


class CompositeLocant:
    """A composite locant of an assembly of three or more rings: ring number and position, cited with the position raised (1⁴)."""

    def __init__(self, ring, position):
        self.ring, self.position = ring, position

    def _key(self):
        return self.ring, float(self.position)

    def __eq__(self, other):
        return isinstance(other, CompositeLocant) and self._key() == other._key()

    def __lt__(self, other):
        return self._key() < other._key()

    def __hash__(self):
        return hash(self._key())

    def __str__(self):
        return superscript_locant(self.ring, self.position)
