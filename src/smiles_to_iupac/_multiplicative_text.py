"""Text helpers for multiplicative names (P-16.5 enclosing marks, P-16.3
multiplying prefixes) shared by the unit, linker and assembly modules.
"""

from ._numerals import numerical_term

_OPENERS = {"(": 0, "[": 1, "{": 2}
_MARKS = ["()", "[]", "{}"]


def enclose(text):
    levels = [_OPENERS[ch] for ch in text if ch in _OPENERS]
    level = (max(levels) + 1) % 3 if levels else 0
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
