"""Ionic skeletal replacement prefixes (P-15.4.1.4, P-72.4 method 2, P-73.4 method 2): a charged skeletal heteroatom is
cited with 'onia', 'ylia', 'ida' or 'uida' in the heteroatom prefixes instead of the neutral 'a' prefix plus the
ium, ylium, ide or uide suffix of the name that gives the PIN. The prefix follows from the bonding number of the atom."""

import re
from collections import Counter

from rdkit import Chem

from ._anion_center import _PARENTS
from ._common import UnsupportedStructure, adjacency
from ._diester_ring_diyl import _system_of, evaluate_skeleton
from ._numerals import alkane_name

_NEUTRAL = {
    5: "bora", 7: "aza", 8: "oxa", 9: "fluora", 13: "alumina", 14: "sila", 15: "phospha", 16: "thia", 17: "chlora",
    31: "galla", 32: "germa", 33: "arsa", 34: "selena", 35: "broma", 49: "inda", 50: "stanna", 51: "stiba", 52: "tellura",
    53: "ioda", 81: "thalla", 82: "plumba", 83: "bisma",
}
_ELEMENT = {prefix: z for z, prefix in _NEUTRAL.items()}
_CHARGE = {"ium": 1, "ylium": 1, "ide": -1, "uide": -1}
_HYDRON = {"ium": 1, "uide": 1, "ylium": -1, "ide": -1}
_COUNT = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
_PREFIX = "|".join(sorted(_ELEMENT, key=len, reverse=True))
_STEMS = "|".join(sorted({alkane_name(k)[:-3] for k in range(1, 100)}, key=len, reverse=True))
_SUFFIX = re.compile(r"^(?P<head>.+?)-(?P<loc>\d+(?:,\d+)*)-(?P<mult>di|tri)?(?P<kind>ium|ylium|uide|ide)$")
_ITEM = re.compile(
    rf"(?P<loc>\d+(?:λ\d+)?(?:,\d+(?:λ\d+)?)*)-(?P<mult>di|tri|tetra|penta|hexa)?(?P<pre>{_PREFIX})"
    rf"(?=-\d|[a-z]*cyclo(?:\[|(?:{_STEMS})an)|[a-z]*spiro\[|(?:{_STEMS})ane?$)"
)
_LOCANT = re.compile(r"(\d+)(?:λ(\d+))?")

# (suffix, group, index of the skeleton atom, indices of the atoms the suffix stands for), most senior first (P-41)
_CLASSES = tuple(
    (suffix, Chem.MolFromSmarts(smarts), anchor, replaced)
    for suffix, smarts, anchor, replaced in (
        ("carboxylic acid", "[#6][CX3](=O)[OX2H1]", 0, (1, 2, 3)),
        ("carbonitrile", "[#6][CX2]#[NX1]", 0, (1, 2)),
        ("carbaldehyde", "[#6][CX3H1]=O", 0, (1, 2)),
        ("one", "[#6][CX3](=O)[#6]", 1, (2,)),
        ("ol", "[#6][OX2H1]", 0, (1,)),
        ("thiol", "[#6][SX2H1]", 0, (1,)),
        ("amine", "[#6][NX3H2]", 0, (1,)),
    )
)
_UNSUPPORTED_SENIOR = tuple(
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=[O,S,N])[N,O,S,F,Cl,Br,I;!$([OX2H1])]",
        "[SX4,SX3,PX4,PX3,AsX4,B](=O)~[O,N]",
    )
)
_AMINE_NITROGEN = Chem.MolFromSmarts("[#6][NX3;+0;!$(N-[!#6;!#1]);!$(N-C=[O,S,N])]")


def ionic_prefix_names(mol, name_parent):
    """[(name, rule)] naming the charged skeletal heteroatoms of the PIN by ionic prefixes, or [] when the PIN holds no
    ionic suffix on a skeleton named by replacement."""
    if not any(a.GetFormalCharge() for a in mol.GetAtoms()):
        return []
    words = name_parent(Chem.MolToSmiles(mol)).split(" ")
    for index, word in enumerate(words):
        pieces = word.split("—")
        for position, piece in enumerate(pieces):
            suffix = _SUFFIX.match(piece)
            if suffix is None:
                continue
            name, charge = _name(mol, suffix)
            pieces[position] = name
            words[index] = "—".join(pieces)
            return [(" ".join(words), "P-73.4" if charge > 0 else "P-72.4")]
    return []


def _name(mol, suffix):
    charge = _CHARGE[suffix.group("kind")]
    try:
        return _from_structure(mol, charge), charge
    except UnsupportedStructure as structured:
        try:
            return _from_name(mol, suffix)
        except UnsupportedStructure:
            raise structured from None


def _from_name(mol, suffix):
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

    def describe(locant, z, lam):
        if locant in charged:
            return z, (lam or _PARENTS[z][1]) + _HYDRON[kind], _CHARGE[kind]
        return None

    text, seen, atoms = _ionize(head, describe)
    if seen != charged:
        raise UnsupportedStructure("a charged locant has no replacement prefix")
    if not _matches(mol, atoms):
        raise UnsupportedStructure("the charged atoms of the structure differ from those the name describes")
    return text, _CHARGE[kind]


def _ionize(head, describe):
    """`head` with the replacement prefix of every charged skeletal atom turned into its ionic prefix; `describe(locant,
    element, lambda)` is (element, bonding number, charge) of a charged atom and None for a neutral one. A prefix group
    that holds charged and neutral atoms of one element is split, the ionic prefixes first."""
    seen, atoms, pieces, last = set(), [], [], 0
    for item in _ITEM.finditer(head):
        z = _ELEMENT[item.group("pre")]
        charged, neutral = {}, []
        for entry in item.group("loc").split(","):
            parsed = _LOCANT.fullmatch(entry)
            locant = int(parsed.group(1))
            described = describe(locant, z, int(parsed.group(2)) if parsed.group(2) else None)
            if described is None:
                neutral.append(entry)
                continue
            if described[0] != z:
                raise UnsupportedStructure("the element of a charged atom differs from its replacement prefix")
            text, lam = _ionic(*described)
            charged.setdefault((text, lam), []).append(locant)
            seen.add(locant)
            atoms.append(described)
        if not charged:
            continue
        parts = []
        for (text, lam), locants in sorted(charged.items(), key=lambda kv: min(kv[1])):
            marked = ",".join(f"{x}λ{lam}" if lam else str(x) for x in sorted(locants))
            parts.append(f"{marked}-{_COUNT[len(locants)]}{text}")
        if neutral:
            parts.append(f"{','.join(neutral)}-{_COUNT[len(neutral)]}{item.group('pre')}")
        pieces += [head[last:item.start()], "-".join(parts)]
        last = item.end()
    return "".join(pieces) + head[last:], seen, atoms


def _from_structure(mol, charge):
    """The ionic-prefix name built from the ring system that holds the charged heteroatoms, numbered as the neutral
    skeleton is, with the senior characteristic group as the suffix (P-73.4: the prefixes are used as neutral ones)."""
    fragment, charged = _cation_fragment(mol, charge)
    found = _system_of(fragment, charged[0])
    if found is None or not set(charged) <= found[1]:
        raise UnsupportedStructure("the charged atoms are not all skeletal atoms of one ring system")
    rings, system = found
    suffix, attach, leaving = _principal_group(fragment, system)
    named = evaluate_skeleton(fragment, adjacency(fragment), "ring", rings, system, attach, leaving, suffix)
    if named is None:
        raise UnsupportedStructure("this ring system has no supported numbering")
    described = {named[2][a]: _describe(fragment, a) for a in charged}
    text, seen, _ = _ionize(named[1], lambda locant, z, lam: described.get(locant))
    if seen != set(described):
        raise UnsupportedStructure("a charged atom has no replacement prefix in the name of its ring system")
    return text


def _cation_fragment(mol, charge):
    holders = [
        atoms for atoms in Chem.GetMolFrags(mol)
        if any(mol.GetAtomWithIdx(a).GetFormalCharge() * charge > 0 and mol.GetAtomWithIdx(a).IsInRing() for a in atoms)
    ]
    if len(holders) != 1:
        raise UnsupportedStructure("the charged skeleton is not in exactly one component")
    editable = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - set(holders[0]), reverse=True):
        editable.RemoveAtom(index)
    fragment = editable.GetMol()
    Chem.SanitizeMol(fragment)
    charged = [a.GetIdx() for a in fragment.GetAtoms() if a.GetFormalCharge()]
    if any(fragment.GetAtomWithIdx(a).GetAtomicNum() not in _NEUTRAL for a in charged) or any(
        fragment.GetAtomWithIdx(a).GetFormalCharge() * charge < 0 for a in charged
    ):
        raise UnsupportedStructure("a charged atom has no replacement prefix or the charges differ in sign")
    return fragment, charged


def _describe(fragment, index):
    kekule = Chem.Mol(fragment)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Chem.KekulizeException:
        pass
    atom = kekule.GetAtomWithIdx(index)
    bonding = int(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs())
    return atom.GetAtomicNum(), bonding, atom.GetFormalCharge()


def _principal_group(fragment, skeleton):
    """(suffix, skeleton atoms bearing it, atoms the suffix stands for) of the senior class present, which must stand on
    the skeleton for the ring system to be the parent; ('', [], set()) when no class outranks the prefixes."""
    if any(fragment.HasSubstructMatch(pattern) for pattern in _UNSUPPORTED_SENIOR):
        raise UnsupportedStructure("a characteristic group senior to the alcohols that the suffix naming does not cover")
    for suffix, pattern, anchor, replaced in _CLASSES:
        matches = fragment.GetSubstructMatches(pattern)
        if suffix == "amine":
            amines = {m[1] for m in fragment.GetSubstructMatches(_AMINE_NITROGEN) if m[1] not in skeleton}
            if amines - {m[1] for m in matches}:
                raise UnsupportedStructure("an amine that is not a primary amino group on the skeleton")
        if not matches:
            continue
        if any(m[anchor] not in skeleton for m in matches):
            raise UnsupportedStructure("the principal characteristic group is not on the charged ring system")
        return suffix, sorted({m[anchor] for m in matches}), {m[i] for m in matches for i in replaced}
    return "", [], set()


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


def _ionic(z, bonding, charge):
    """(prefix, lambda number to keep) of a charged atom of element `z` with the given bonding number: the cationic
    prefix of a bonding number one above ('onia') or below ('ylia') the standard one, the anionic one one below
    ('ida') or above ('uida') it."""
    stem, standard = _PARENTS[z]
    if charge > 0:
        if bonding == standard + 1:
            return ("bismuthonia" if z == 83 else _NEUTRAL[z][:-1] + "onia"), None
        if bonding == standard - 1:
            return f"{stem}ylia", None
        raise UnsupportedStructure("the cationic prefixes cover a bonding number one above or below the standard one")
    if bonding == standard - 1:
        return f"{stem}ida", None
    if bonding == standard + 1:
        return f"{stem}uida", None
    if bonding > standard + 1:
        return f"{stem}uida", bonding - 1
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
