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
    r"|\[\d+(?:,\d+)*-[a-z]{1,2}['\u2032]*(?::\d+(?:,\d+)*-[a-z]{1,2}['\u2032]*)+\]"  # multiparent: [1,2-b:5,4-b']
)


_ADDED_HYDROGEN = re.compile(r"\(\d+[a-z]?H\)")


def nesting_text(text):
    """`text` without the marks that P-16.5.4.1 leaves out of the nesting order: the descriptor brackets of ring
    assembly, von Baeyer, spiro and fusion names, and the parentheses of added hydrogen."""
    return _ADDED_HYDROGEN.sub("", _ASSEMBLY_BRACKETS.sub("", text))


_CODE = r"(?:\d+[a-z]?['\u2032]*)?(?:[RSEZrsMPmp]|seqCis|seqTrans)"
_DESCRIPTOR = re.compile(rf"\((?:{_CODE}(?:,{_CODE})*|(?:T|SP|SS|TBPY|OC|SPY)-[\dA-Z\u2032-]+)\)-")


def _mark_index(core):
    """Index into ( [ { of the marks that enclose `core`: the one above its deepest group. A stereodescriptor opening a
    group is no group of its own, but a group whose mark would repeat its first character takes the next one
    (P-16.5.4)."""
    stack, deepest, i = [], 0, 0
    lead = _DESCRIPTOR.match(core)
    if lead:
        i = lead.end()
    while i < len(core):
        ch = core[i]
        if ch in _OPENERS:
            stack.append((i, 0))
            lead = _DESCRIPTOR.match(core, i + 1)
            i = lead.end() if lead else i + 1
            continue
        if ch in ")]}" and stack:
            opened, below = stack.pop()
            rank = below + 1
            while core[opened + 1:opened + 2] == _MARKS[(rank - 1) % 3][0]:
                rank += 1
            if stack:
                stack[-1] = (stack[-1][0], max(stack[-1][1], rank))
            else:
                deepest = max(deepest, rank)
        i += 1
    return deepest % 3


def enclose(text):
    core = nesting_text(text)
    level = _mark_index(core)
    if core.startswith(_MARKS[level][0]):
        level = (level + 1) % 3
    left, right = _MARKS[level]
    return f"{left}{text}{right}"


def join_parts(pieces):
    """The pieces of a concatenated substituent, with a hyphen between a name and the locants that follow it (P-16.2.4.1)."""
    text = ""
    for piece in pieces:
        text += ("-" if text and text[-1] not in "([{-" and piece[:1].isdigit() else "") + piece
    return text


def multiplier_word(count, use_bis):
    if count == 1:
        return ""
    if use_bis:
        return {2: "bis", 3: "tris"}.get(count, numerical_term(count) + "kis")
    return numerical_term(count)


def primed_locants(locant, count):
    return ",".join(str(locant) + "'" * i for i in range(count))


_NUMERIC_STEMS = (
    "pentan", "hexan", "heptan", "octan", "nonan",
    "dec", "undec", "dodec", "tridec", "tetradec", "pentadec", "hexadec", "heptadec", "octadec", "nonadec", "icos",
)


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


class JunctionBondLocant(PrimedLocant):
    """Locant of a double bond joining two components of an assembly, cited as 1(1′) (P-93.5.1.4.2.2)."""

    def __init__(self, first, second):
        super().__init__(first.primes, first.number)
        self.second = second

    def __str__(self):
        return f"{PrimedLocant.__str__(self)}({self.second})"


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
