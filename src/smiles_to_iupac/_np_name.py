"""Names built on a stereoparent candidate: modified parent, groups, unsaturation and configuration (P-101.3 - P-101.7)."""

from collections import Counter

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, halogen_substituents, multiplied_word
from ._np_core import FACES, center_signature, exo_faces, is_numbered, loc_key, orientation_sign, parent_h
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_SUFFIX = {"ketone": "one", "alcohol": "ol", "amine": "amine", "ester_o": "yl", "aldehyde": "al", "diyl": "diyl", "yl": "yl"}
_ACYL_SUFFIX = {
    ("acid", "c"): "carboxylic acid",
    ("acid", "o"): "oic acid",
    ("ester", "c"): "carboxylate",
    ("ester", "o"): "oate",
    ("amide", "c"): "carboxamide",
    ("amide", "o"): "amide",
}
_PREFIX = {"alcohol": "hydroxy", "ketone": "oxo", "amine": "amino", "aldehyde": "oxo"}
_SENIORITY = ["acid", "ester", "ester_o", "diyl", "amide", "aldehyde", "ketone", "alcohol", "amine"]
_COUNTED = {
    "acid": "[CX3](=O)[OX2H1]",
    "amide": "[CX3;!R](=O)[NX3;!R]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "alcohol": "[OX2H][#6]",
}
_ACYL_CLASSES = ("acid", "ester", "ester_o", "amide", "diyl")
_LETTERS = "abcdefghij"


class _IntLoc(int):
    def __new__(cls, locant, face=""):
        text = str(locant)
        digits = text.rstrip("′″‴")
        primes = len(text) - len(digits)
        obj = super().__new__(cls, int(digits) + 1000 * primes)
        obj.base = text
        obj.text = f"{text}{FACES.get(face, '')}"
        return obj

    def __str__(self):
        return self.text

    __repr__ = __str__
    __format__ = lambda self, spec: self.text


class _StrLoc(str):
    def __new__(cls, locant, face=""):
        obj = super().__new__(cls, f"{locant}{FACES.get(face, '')}")
        obj.base = locant
        return obj


def Loc(locant, face=""):
    """A locant carrying a configuration symbol (17β); numeric locants stay ordered as numbers."""
    return _IntLoc(locant, face) if str(locant).rstrip("′″‴").isdigit() else _StrLoc(locant, face)


def final_labels(skel):
    """{provisional label: final locant}: 'homo' atoms (P-101.3.2.2) and unnumbered heteroatoms (P-101.4.3)."""
    homo = [op for op in skel.ops if op[0] == "homo"]
    labels = {}
    letters = Counter()
    keyed = []
    for _, kind, data, label in homo:
        if kind == "terminal":
            base = data
        elif kind == "atomic":
            base = data[0]
        else:
            base = f"{data[0]}({data[1]})"
        keyed.append((loc_key(data if kind == "terminal" else data[0]), base, kind, label))
    for _, base, kind, label in sorted(keyed, key=lambda t: t[:2]):
        key = base if kind != "bond" else base.split("(")[0]
        labels[label] = f"{base}{_LETTERS[letters[key]]}"
        letters[key] += 1
    unnumbered = [a for a in skel.adj if not a.startswith("h") and not is_numbered(a)]
    taken = Counter()
    for atom in sorted(unnumbered, key=loc_key):
        neighbors = sorted((n for n in skel.adj[atom] if is_numbered(n)), key=loc_key)
        if not neighbors:
            continue
        base = neighbors[0]
        labels[atom] = f"{base}{_LETTERS[taken[base]]}"
        taken[base] += 1
    return labels


def _final(skel):
    labels = final_labels(skel)
    return lambda loc: labels.get(loc, loc)


def _carboxyl(view, atom, parent, mapped, ignored=frozenset()):
    others = [n for n in view.adj[atom] if n != parent and n not in ignored]
    oxo = [n for n in others if view.elem[n] == "O" and view.order[frozenset((atom, n))] == 2]
    rest = [n for n in others if n not in oxo]
    if len(oxo) != 1 or len(rest) != 1 or view.order[frozenset((atom, rest[0]))] != 1:
        return None
    x = rest[0]
    mol = view.mol
    if view.elem[x] == "O" and len(view.adj[x]) == 1 and mol.GetAtomWithIdx(x).GetTotalNumHs() == 1:
        return "acid", None
    if view.elem[x] == "O" and len(view.adj[x]) == 2:
        alkyl = next(n for n in view.adj[x] if n != atom)
        if view.elem[alkyl] != "C" or alkyl in mapped:
            raise UnsupportedStructure("a lactone or a non-carbon ester group on a natural product is not supported")
        return "ester", (x, alkyl)
    if view.elem[x] == "N" and len(view.adj[x]) == 1 and mol.GetAtomWithIdx(x).GetTotalNumHs() == 2:
        return "amide", None
    return None


def _acyl_oxygen(view, oxygen, parent, mapped):
    acyl = next((n for n in view.adj[oxygen] if n != parent), None)
    if acyl is None or view.elem[acyl] != "C" or acyl in mapped:
        return None
    oxo = [n for n in view.adj[acyl] if view.elem[n] == "O" and len(view.adj[n]) == 1 and view.order[frozenset((acyl, n))] == 2]
    if len(oxo) != 1 or len(view.adj[acyl]) > 3:
        return None
    if any(n not in (oxygen, oxo[0]) and view.elem[n] != "C" for n in view.adj[acyl]):
        return None
    return acyl


def _acid_anion(view, oxygen, acyl, mapped):
    from .core import smiles_to_iupac

    keep, stack = {acyl}, [acyl]
    while stack:
        for n in view.adj[stack.pop()]:
            if n != oxygen and n not in keep:
                if n in mapped:
                    raise UnsupportedStructure("a ring-closing O-acyl group is not supported")
                keep.add(n)
                stack.append(n)
    editable = Chem.RWMol(view.mol)
    editable.GetAtomWithIdx(oxygen).SetNumExplicitHs(1)
    editable.GetAtomWithIdx(oxygen).SetNoImplicit(True)
    for idx in sorted((i for i in view.adj if i not in keep | {oxygen}), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    name = smiles_to_iupac(Chem.MolToSmiles(acid))
    if not name.endswith("ic acid") or " " in name[: -len(" acid")]:
        raise UnsupportedStructure("the acid part of this ester is not a plain acid")
    return name[: -len("ic acid")] + "ate"


class Groups:
    def __init__(self):
        self.classes = {}
        self.branches = []
        self.rings = []

    def add(self, cls, loc, atoms, extra=None):
        self.classes.setdefault(cls, []).append((loc, atoms, extra or {}))


def classify(cand, view, ignored=frozenset()):
    """Groups on the skeleton atoms of a candidate; raises UnsupportedStructure for shapes not handled."""
    skel = cand.skel
    mapping = cand.mapping
    mapped = set(mapping.values())
    image = {a: loc for loc, a in mapping.items()}
    ring_atoms = skel.ring_atoms()
    groups = Groups()
    mol = view.mol
    for loc, atom in mapping.items():
        outside = [n for n in view.adj[atom] if n not in mapped and n not in ignored]
        leaf = loc not in ring_atoms and len(skel.adj[loc]) == 1
        if leaf and view.elem[atom] == "C":
            parent_atom = mapping[next(iter(skel.adj[loc]))]
            carboxyl = _carboxyl(view, atom, parent_atom, mapped, ignored)
            if carboxyl:
                kind, ester = carboxyl
                groups.add(kind, loc, {atom}, {"kind": "o", "ester": ester})
                continue
        oxo = [n for n in outside if view.elem[n] == "O" and len(view.adj[n]) == 1 and view.order[frozenset((atom, n))] == 2]
        hydroxyl = [n for n in outside if view.elem[n] == "O" and len(view.adj[n]) == 1 and mol.GetAtomWithIdx(n).GetTotalNumHs() == 1]
        amino = [n for n in outside if view.elem[n] == "N" and len(view.adj[n]) == 1 and mol.GetAtomWithIdx(n).GetTotalNumHs() == 2]
        if any(view.order[frozenset((atom, n))] != 1 for n in outside if n not in oxo):
            raise UnsupportedStructure("an exocyclic multiple bond on a natural-product skeleton is not supported")
        if oxo:
            if len(outside) != 1:
                raise UnsupportedStructure("an acyl group on a natural-product skeleton is not supported")
            if loc not in ring_atoms and mol.GetAtomWithIdx(atom).GetTotalNumHs() > 0:
                groups.add("aldehyde", loc, {oxo[0]})
            else:
                groups.add("ketone", loc, {oxo[0]})
            continue
        for n in outside:
            if n in hydroxyl:
                groups.add("alcohol", loc, {n})
            elif n in amino:
                groups.add("amine", loc, {n})
            else:
                carboxyl = _carboxyl(view, n, atom, mapped, ignored) if view.elem[n] == "C" else None
                if carboxyl and loc in ring_atoms or (carboxyl and loc not in ring_atoms):
                    kind, ester = carboxyl
                    groups.add(kind, loc, {n}, {"kind": "c", "ester": ester, "anchor": n, "root": n})
                    continue
                acyl = _acyl_oxygen(view, n, atom, mapped) if view.elem[n] == "O" and len(view.adj[n]) == 2 else None
                if acyl is not None:
                    groups.add("ester_o", loc, {n}, {"anchor": n, "root": n, "acyl": acyl})
                    continue
                groups.branches.append((loc, n))
    _pair_acetals(view, mapping, groups)
    return groups


def _pair_acetals(view, mapping, groups):
    """A carbon bonded to two oxygens that both sit on skeleton atoms closes a cyclic acetal, ketal or carbonate."""
    by_carbon = {}
    for loc, root in groups.branches:
        if view.elem[root] == "O" and len(view.adj[root]) == 2:
            carbon = next(n for n in view.adj[root] if n != mapping[loc])
            by_carbon.setdefault(carbon, []).append((loc, root))
    for carbon, members in by_carbon.items():
        if len(members) != 2 or view.elem[carbon] != "C":
            continue
        others = [n for n in view.adj[carbon] if n not in {r for _, r in members}]
        groups.branches = [b for b in groups.branches if b not in members]
        groups.add("diyl", members[0][0], {r for _, r in members}, {"carbon": carbon, "others": others, "second": members[1][0], "anchors": [r for _, r in members]})


def alkyl_count(cand, view, groups):
    """Number of carbon-chain substituent groups (modifications that compete with 'nor' and 'homo', P-101.3.7.1)."""
    count = 0
    ring_atoms = cand.skel.ring_atoms()
    for loc, root in groups.branches:
        if view.elem[root] != "C":
            continue
        if any(view.elem[n] in "ONS" and view.order[frozenset((root, n))] >= 2 for n in view.adj[root]):
            continue
        if any(view.elem[n] == "N" and view.order[frozenset((root, n))] == 3 for n in view.adj[root]):
            continue
        if loc not in ring_atoms:
            return None
        count += 1
    return count


def principal_groups_left_outside(view, classes, principal):
    """True when a group of the principal class or a senior one sits outside the parent (P-44.1.1.2, maximum number)."""
    senior = _SENIORITY if principal is None else _SENIORITY[: _SENIORITY.index(principal) + 1]
    for cls in senior:
        smarts = _COUNTED.get(cls)
        if smarts is None:
            continue
        total = len(view.mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
        if total > len(classes.get(cls, [])):
            return True
    return False
