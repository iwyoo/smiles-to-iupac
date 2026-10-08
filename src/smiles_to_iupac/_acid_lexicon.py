"""Names and seniority of the acid groups of P-65 and their functional-replacement analogues (P-43.1, Tables 4.2 and 4.3).
An acid group is a centre (C, S, Se, Te) with =X positions and a -Y-H chain; replaced oxygens add infixes cited
alphabetically, the last one carrying 'ic acid', with italic letters for the tautomer when several chalcogens occur.
"""

import re
from dataclasses import dataclass

CHALCOGENS = ("O", "S", "Se", "Te")
_REPLACED = (*CHALCOGENS, "NH", "NNH2")
_INFIX = {"S": "thio", "Se": "seleno", "Te": "telluro", "NH": "imido", "NNH2": "hydrazono"}
_MULTIPLIER = {1: "", 2: "di", 3: "tri", 4: "tetra"}
_PLAIN_WORD = {
    ("C", 1): "carboxylic acid",
    ("S", 2): "sulfonic acid",
    ("S", 1): "sulfinic acid",
    ("Se", 2): "selenonic acid",
    ("Se", 1): "seleninic acid",
    ("Te", 2): "telluronic acid",
    ("Te", 1): "tellurinic acid",
}
_STEM = {("S", 2): "sulfon", ("S", 1): "sulfin", ("Se", 2): "selenon", ("Se", 1): "selenin", ("Te", 2): "telluron", ("Te", 1): "tellurin"}
_CENTER_RANK = {("C", 1): 0, ("S", 2): 1, ("S", 1): 2, ("Se", 2): 3, ("Se", 1): 4, ("Te", 2): 5, ("Te", 1): 6}
_NITROGEN_FAMILY = {"NH": 1, "NNH2": 2}


@dataclass(frozen=True)
class AcidSpec:
    center: str
    oxo: tuple
    y: tuple
    anion: bool = False

    @property
    def kind(self):
        return (self.center, len(self.oxo))

    @property
    def plain(self):
        return not self.anion and self.oxo == ("O",) * len(self.oxo) and self.y == ("O",)

    @property
    def key(self):
        return f"acid:{self.center}:{'+'.join(self.oxo)}:{'.'.join(self.y)}{'!' if self.anion else ''}"


def make_spec(center, oxo, y, anion=False):
    return AcidSpec(center, tuple(sorted(oxo, key=_REPLACED.index)), tuple(y), anion)


def spec_from_key(key):
    _, center, oxo, y = key.split(":")
    anion = y.endswith("!")
    return AcidSpec(center, tuple(oxo.split("+")), tuple(y.rstrip("!").split(".")), anion)


def rank_key(spec):
    """Seniority of an acid group, smaller is senior (P-41 class 7a, P-43.1)."""
    nitrogen = max((_NITROGEN_FAMILY.get(x, 0) for x in spec.oxo), default=0)
    unmodified = all(x == "O" or x in _NITROGEN_FAMILY for x in spec.oxo) and spec.y == ("O",)
    atoms = (*spec.oxo, *spec.y)
    return (
        0 if spec.anion else 1,
        _CENTER_RANK[spec.kind],
        nitrogen,
        0 if unmodified else 1,
        -len(spec.y),
        tuple(-atoms.count(e) for e in _REPLACED),
        tuple(-spec.y.count(e) for e in CHALCOGENS),
        CHALCOGENS.index(spec.y[-1]),
    )


def _peroxo_word(y):
    replaced = {}
    for atom in y:
        if atom != "O":
            replaced[_INFIX[atom]] = replaced.get(_INFIX[atom], 0) + 1
    return "".join(_MULTIPLIER[replaced[w]] + w for w in sorted(replaced)) + "peroxo"


def _terms(spec):
    """Infixes of an acid group as sorted (word, count, compound) triples."""
    simple = {}
    for atom in (*spec.oxo, *(spec.y if len(spec.y) == 1 else ())):
        if atom != "O":
            simple[_INFIX[atom]] = simple.get(_INFIX[atom], 0) + 1
    terms = [(word, count, False) for word, count in simple.items()]
    if len(spec.y) == 2:
        word = _peroxo_word(spec.y)
        terms.append((word, 1, word != "peroxo"))
    return sorted(terms)


def _term_text(term, ending):
    word, count, compound = term
    text = _MULTIPLIER[count] + word
    if ending == "ic":
        if word in ("imido", "hydrazono"):
            text = text[:-1]
        text += "ic"
    elif ending == "yl":
        text += "yl"
    return f"({text})" if compound else text


def _infixes(spec, ending):
    terms = _terms(spec)
    texts = []
    for i, term in enumerate(terms):
        text = _term_text(term, ending if i == len(terms) - 1 else None)
        # P-16.7.1: the terminal 'o' of an infix is elided before the vowel of the next one
        if i < len(terms) - 1 and text.endswith("o") and terms[i + 1][0][0] in "aeiou":
            text = text[:-1]
        texts.append(text)
    return "".join(texts)


def tautomer_letters(spec):
    """Element letters of the -Y-H chain when several chalcogens are present (P-65.1.5.1)."""
    if len(spec.y) == 2:
        return "".join(spec.y) if spec.y[0] != spec.y[1] else ""
    elements = {x for x in spec.oxo if x in CHALCOGENS} | set(spec.y)
    return "".join(spec.y) if len(elements) > 1 else ""


def _acid_tail(spec, count=1):
    letters = tautomer_letters(spec)
    return f" {','.join([letters] * count)}-acid" if letters else " acid"


_ANION_ENDING = re.compile(r"ic(?: [A-Za-z,]+-acid| acid)?(\)?)$")


_LETTERED_TAIL = re.compile(r"\((\w+?)ic\) ([A-Za-z,]+)-acid$")


def _anionic(text):
    lettered = _LETTERED_TAIL.search(text)
    if lettered:
        return text[: lettered.start()] + f"({lettered.group(2)}-{lettered.group(1)}ate)"
    return _ANION_ENDING.sub(lambda m: "ate" + (m.group(1) or ""), text)


def chain_suffix(spec, count):
    """(body, tail) of the suffix of an acid group whose carbon is part of an
    acyclic chain ('oic acid' type): the body fuses with the chain stem."""
    if spec.anion:
        body, tail = chain_suffix(AcidSpec(spec.center, spec.oxo, spec.y), count)
        return _anionic(body + tail), ""
    if spec.plain:
        return {1: "oic", 2: "dioic", 3: "trioic", 4: "tetraoic"}[count], " acid"
    text = _infixes(spec, "ic")
    tail = _acid_tail(spec, count)
    first = _terms(spec)[0]
    if count > 1:
        # 'dithioic' would read as one group with two sulfur atoms, so thio-type infixes need 'bis'
        if first[0] in ("thio", "seleno", "telluro") or first[2]:
            return f"{'bis' if count == 2 else 'tris'}({text}{tail})", ""
        return _MULTIPLIER[count] + text, tail
    if first[1] > 1 and tail == " acid":
        return f"({text} acid)", ""
    return text, tail


def _head(spec, text):
    if spec.center == "C":
        return "carbox" if text[0] in "aeiouy" else "carbo"
    stem = _STEM[spec.kind]
    return stem if text[0] in "aeiouy" else stem + "o"


def carbo_suffix(spec, count):
    """Suffix of an acid group on a ring or heteroatom parent ('carboxylic acid' type)."""
    if spec.anion:
        return _anionic(carbo_suffix(AcidSpec(spec.center, spec.oxo, spec.y), count))
    if spec.plain:
        return _MULTIPLIER[count] + _PLAIN_WORD[spec.kind]
    text = _infixes(spec, "ic")
    return _MULTIPLIER[count] + _head(spec, text) + text + _acid_tail(spec)


def acyl_suffix(spec, chain, count):
    """Ending of the acyl group of an acid (-Y-H removed), cited after the parent name."""
    oxo_only = make_spec(spec.center, spec.oxo, ("O",))
    terms = _terms(oxo_only)
    if not terms:
        if chain:
            return {1: "oyl", 2: "dioyl", 3: "trioyl", 4: "tetraoyl"}[count]
        stem = "carbon" if spec.center == "C" else _STEM[spec.kind]
        return _MULTIPLIER[count] + stem + "yl"
    text = "".join(_term_text(t, "yl" if i == len(terms) - 1 else None) for i, t in enumerate(terms))
    if chain:
        return text if count == 1 else f"bis({text})"
    return _MULTIPLIER[count] + _head(oxo_only, text) + text
