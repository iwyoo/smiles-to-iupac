"""Ionic skeletal replacement prefixes (P-15.4.1.4, P-72.4 method 2, P-73.4 method 2): a charged skeletal heteroatom is
cited with 'onia', 'ylia', 'ida' or 'uida' in the heteroatom prefixes instead of the neutral 'a' prefix plus the
ium, ylium, ide or uide suffix of the name that gives the PIN. The prefix follows from the bonding number of the atom."""

import re
from collections import Counter

from rdkit import Chem

from ._anion_center import _PARENTS
from ._common import UnsupportedStructure
from ._numerals import alkane_name

_NEUTRAL = {
    5: "bora", 7: "aza", 8: "oxa", 9: "fluora", 13: "alumina", 14: "sila", 15: "phospha", 16: "thia", 17: "chlora",
    31: "galla", 32: "germa", 33: "arsa", 34: "selena", 35: "broma", 49: "inda", 50: "stanna", 51: "stiba", 52: "tellura",
    53: "ioda", 81: "thalla", 82: "plumba", 83: "bisma",
}
_ELEMENT = {prefix: z for z, prefix in _NEUTRAL.items()}
_CHARGE = {"ium": 1, "ylium": 1, "ide": -1, "uide": -1}
_HYDRON = {"ium": 1, "uide": 1, "ylium": -1, "ide": -1}
_PREFIX = "|".join(sorted(_ELEMENT, key=len, reverse=True))
_STEMS = "|".join(sorted({alkane_name(k)[:-3] for k in range(1, 100)}, key=len, reverse=True))
_SUFFIX = re.compile(r"^(?P<head>.+?)-(?P<loc>\d+(?:,\d+)*)-(?P<mult>di|tri)?(?P<kind>ium|ylium|uide|ide)$")
_ITEM = re.compile(
    rf"(?P<loc>\d+(?:λ\d+)?(?:,\d+(?:λ\d+)?)*)-(?P<mult>di|tri|tetra|penta|hexa)?(?P<pre>{_PREFIX})"
    rf"(?=-\d|[a-z]*cyclo(?:\[|(?:{_STEMS})an)|[a-z]*spiro\[|(?:{_STEMS})ane?$)"
)


def ionic_prefix_names(mol, name_parent):
    """[(name, rule)] naming the charged skeletal heteroatoms of the PIN by ionic prefixes, or [] when the PIN holds no
    ionic suffix on a skeleton named by replacement."""
    if not any(a.GetFormalCharge() for a in mol.GetAtoms()):
        return []
    words = name_parent(Chem.MolToSmiles(mol)).split(" ")
    for index, word in enumerate(words):
        suffix = _SUFFIX.match(word)
        if suffix is None:
            continue
        name, charge = _rewrite(mol, suffix)
        return [(" ".join(words[:index] + [name] + words[index + 1:]), "P-73.4" if charge > 0 else "P-72.4")]
    return []


def _rewrite(mol, suffix):
    kind, head = suffix.group("kind"), suffix.group("head")
    charged = {int(x) for x in suffix.group("loc").split(",")}
    if suffix.group("mult") is None:
        if not head.endswith("an"):
            raise UnsupportedStructure("the ionic suffix does not follow an alkane or cycloalkane ending")
        head += "e"
    items = list(_ITEM.finditer(head))
    if items:
        hetero = {int(x.split("λ")[0]) for m in items for x in m.group("loc").split(",")}
        _require_plain_prefixes(head[: items[0].start()], hetero)
    seen, atoms, pieces, last = set(), [], [], 0
    for item in items:
        locants = [re.fullmatch(r"(\d+)(?:λ(\d+))?", x) for x in item.group("loc").split(",")]
        here = {int(m.group(1)) for m in locants}
        if not here & charged:
            continue
        if not here <= charged:
            raise UnsupportedStructure("charged and uncharged heteroatoms of one element share a prefix group")
        z = _ELEMENT[item.group("pre")]
        results = {_ionic(z, kind, int(m.group(2)) if m.group(2) else _PARENTS[z][1]) for m in locants}
        if len(results) != 1:
            raise UnsupportedStructure("the heteroatoms of one prefix group need different ionic prefixes")
        text, lam, n = results.pop()
        marked = ",".join(f"{m.group(1)}λ{lam}" if lam else m.group(1) for m in locants)
        pieces += [head[last:item.start()], f"{marked}-{item.group('mult') or ''}{text}"]
        last = item.end()
        seen |= here
        atoms += [(z, n, _CHARGE[kind])] * len(here)
    if seen != charged:
        raise UnsupportedStructure("a charged locant has no replacement prefix")
    if not _matches(mol, atoms):
        raise UnsupportedStructure("the charged atoms of the structure differ from those the name describes")
    return "".join(pieces) + head[last:], _CHARGE[kind]


def _require_plain_prefixes(text, hetero_locants):
    """Without the ionic suffix, a hydroxy, amino, carboxy or similar prefix would turn into the principal characteristic
    group: only hydrocarbyl, alkoxy, halo and oxo on a heteroatom stay prefixes."""
    for oxo in re.finditer(r"(\d+(?:,\d+)*)-(?:di|tri|tetra)?oxo", text):
        if not {int(x) for x in oxo.group(1).split(",")} <= hetero_locants:
            raise UnsupportedStructure("an oxo prefix on a carbon atom becomes a suffix without the ionic suffix")
    bare = re.sub(r"\((?:\d+[A-Za-z]?(?:,\d+[A-Za-z]?)*)?[RSEZrsez]+\)", "", text)
    bare = re.sub(r"[\d,\-()\[\]{}λ]|\b[NO]\b", "", bare)
    allowed = rf"(?:di|tri|tetra|penta|hexa|bis|tris|tert|cyclo|phenyl|(?:{_STEMS})(?:an)?yl|(?:{_STEMS})oxy|oxo|fluoro|chloro|bromo|iodo)*"
    if not re.fullmatch(allowed, bare):
        raise UnsupportedStructure("a prefix other than hydrocarbyl, alkoxy, halo or oxo would change the name without the ionic suffix")


def _ionic(z, kind, parent_bonding):
    """(prefix, lambda number to keep, bonding number) of a charged atom formed from the parent hydride with
    `parent_bonding` bonds by adding a hydron or a hydride, or removing a hydride or a hydron."""
    stem, standard = _PARENTS[z]
    n = parent_bonding + _HYDRON[kind]
    if kind in ("ium", "ylium"):
        if n == standard + 1:
            return ("bismuthonia" if z == 83 else _NEUTRAL[z][:-1] + "onia"), None, n
        if n == standard - 1:
            return f"{stem}ylia", None, n
        raise UnsupportedStructure("the cationic prefixes cover a bonding number one above or below the standard one")
    if n == standard - 1:
        return f"{stem}ida", None, n
    if n == standard + 1:
        return f"{stem}uida", None, n
    if n > standard + 1:
        return f"{stem}uida", n - 1, n
    raise UnsupportedStructure("this bonding number is not covered by an anionic prefix")


def _matches(mol, expected):
    """Whether the charged atoms of one component of `mol` are exactly the (element, bonding number, charge) triples."""
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Chem.KekulizeException:
        pass
    wanted = Counter(expected)
    for fragment in Chem.GetMolFrags(kekule, asMols=True, sanitizeFrags=False):
        fragment.UpdatePropertyCache(strict=False)
        found = Counter(
            (a.GetAtomicNum(), int(sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) + a.GetTotalNumHs()), a.GetFormalCharge())
            for a in fragment.GetAtoms()
            if a.GetFormalCharge()
        )
        if found == wanted:
            return True
    return False
