"""Amplificants of phane parent hydrides (P-26.2.2, P-26.5.2): a ring or ring system cited by its parent hydride name
('pyridine' -> 'pyridina', '1,2-oxazole' -> '[1,2]oxazola') and numbered by the rules of that parent. A heteromonocycle
with more than ten members or a heteroatom-bearing von Baeyer or spiro system is cited as the carbon skeleton with
skeletal replacement ('a') heteroatoms in composite locants."""

import re
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import BondType, RWMol

from ._common import UnsupportedStructure, adjacency, superscript_locant
from ._multiplicative_groups import ring_seniority_key
from ._ring_diyl_numbering import _PREFIX, _bare_skeleton, is_hydro_fusion_system, monocycle_numberings, system_numberings
from ._steroid_named import _embeddings as _steroid_embeddings

_INDICATED_HYDROGEN = re.compile(r"^(\d+H(?:,\d+H)*)-")
_LEADING_LOCANTS = re.compile(r"^(\d+[a-z]?(?:,\d+[a-z]?)*)-(?=[a-z])")
_MULTIPLYING_START = ("bi", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca", "di")
_LARGEST_HANTZSCH_WIDMAN_RING = 10
def _build_morphinan_query():
    rw = RWMol()
    idx = {i: rw.AddAtom(Chem.Atom(7 if i == 0 else 6)) for i in range(17)}
    for i in (7, 8, 9, 10, 11, 12):
        rw.GetAtomWithIdx(idx[i]).SetIsAromatic(True)
    single = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (5, 6), (6, 7), (12, 3), (3, 13), (13, 14), (14, 15), (15, 16), (16, 4)]
    ring = [(7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 7)]
    for a, b in single:
        rw.AddBond(idx[a], idx[b], BondType.SINGLE)
    for a, b in ring:
        rw.AddBond(idx[a], idx[b], BondType.AROMATIC)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol, idx


_MORPHINAN, _MORPHINAN_INDEX = _build_morphinan_query()
_MORPHINAN_LOCANTS = {0: 17, 1: 16, 2: 15, 3: 13, 4: 14, 5: 9, 6: 10, 7: 11, 8: 1, 9: 2, 10: 3, 11: 4, 12: 12, 13: 5, 14: 6, 15: 7, 16: 8}


class AmpLoc(int):
    """Locant of an amplificant atom in its own numbering: 1..8 or, for a fusion atom, '4a'."""

    def __new__(cls, text):
        digits, letter, interior = re.fullmatch(r"(\d+)([a-h]?)(\d*)", str(text)).groups()
        order = int(digits) * 1000 + (ord(letter) - 96 if letter else 0) * 100 + (int(interior) if interior else 0)
        obj = super().__new__(cls, order)
        obj.text = str(text)
        return obj

    def __str__(self):
        return self.text

    __repr__ = __str__
    __format__ = lambda self, spec: self.text


class PhaneLoc(int):
    """Locant of a phane skeleton atom or, with `local`, of an amplificant atom (primary 1, local 4 -> '1⁴')."""

    def __new__(cls, primary, local=None):
        obj = super().__new__(cls, primary * 1000000 + (int(local) if local is not None else 0))
        obj.primary, obj.local = primary, local
        return obj

    def __str__(self):
        return superscript_locant(self.primary, self.local) if self.local is not None else str(self.primary)

    __repr__ = __str__
    __format__ = lambda self, spec: str(self)


@dataclass
class Numbered:
    local: dict
    ih: tuple = ()
    hydro: tuple = ()


@dataclass
class Amplificant:
    atoms: frozenset
    prefix: str
    rank: tuple
    mancude: bool
    numberings: list = field(default_factory=list)
    replaced: dict = field(default_factory=dict)
    added: object = None
    bonds: list = field(default_factory=list)
    saturated_parent: bool = False


def needs_bis(prefix):
    return prefix[0] in "[0123456789" or prefix.startswith(_MULTIPLYING_START)


def amplification_prefix(parent):
    """('prefix', ih text) of a parent hydride name: final 'e' becomes 'a' (or 'a' is added), heteroatom locants are
    bracketed (P-26.2.2.1, P-26.0)."""
    stem = _INDICATED_HYDROGEN.sub("", parent)
    leading = _LEADING_LOCANTS.match(stem)
    if leading:
        stem = f"[{leading.group(1)}]" + stem[leading.end():]
    return (stem[:-1] if stem.endswith("e") else stem) + "a"


def _ring_order(graph, atoms):
    start = min(atoms)
    cycle, previous = [start], None
    while len(cycle) < len(atoms):
        following = next(n for n in sorted(graph[cycle[-1]]) if n in atoms and n != previous and n not in cycle)
        previous = cycle[-1]
        cycle.append(following)
    return cycle


def _saturated_analogue(bare, carbon):
    """The bare ring system with every bond single (and every atom carbon when `carbon`): the saturated parent whose
    name an amplification prefix is derived from when 'ene' endings or skeletal replacement carry the rest."""
    editable = Chem.RWMol(bare)
    for atom in editable.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
        if carbon:
            atom.SetAtomicNum(6)
    for bond in editable.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    analogue = editable.GetMol()
    Chem.SanitizeMol(analogue)
    return analogue


def _is_saturated(bare):
    return not any(a.GetIsAromatic() for a in bare.GetAtoms()) and all(
        b.GetBondTypeAsDouble() == 1.0 for b in bare.GetBonds()
    )


def _numberings_of(bare):
    graph = adjacency(bare)
    atoms = set(range(bare.GetNumAtoms()))
    rings = [tuple(r) for r in bare.GetRingInfo().AtomRings()]
    if bare.GetNumBonds() - bare.GetNumAtoms() + 1 == 1:
        return monocycle_numberings(bare, _ring_order(graph, atoms), set(), 0)
    return system_numberings(bare, graph, rings, atoms)


def _parent_name(numbering, bare):
    if bare.GetNumAtoms() == 6 and all(a.GetIsAromatic() and a.GetAtomicNum() == 6 for a in bare.GetAtoms()):
        return "benzene"
    return numbering.parent_stem or numbering.text([], 0, suffix="")


def _added_hydrogen_ring(mol, atoms, free_atom):
    """(aromatic ring analogue, atom with added hydrogen) of a six-membered ring entered by an ylidene bond
    (P-29.3.4.1), else None."""
    ring = sorted(atoms)
    if len(ring) != 6 or free_atom not in atoms:
        return None
    members = [mol.GetAtomWithIdx(a) for a in ring]
    if any(a.GetAtomicNum() not in (6, 7) for a in members) or sum(a.GetAtomicNum() == 7 for a in members) > 1:
        return None
    doubles = sum(
        1 for b in mol.GetBonds() if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms and b.GetBondTypeAsDouble() == 2.0
    )
    sp3 = [
        a.GetIdx()
        for a in members
        if a.GetIdx() != free_atom
        and not any(b.GetBondTypeAsDouble() == 2.0 for b in a.GetBonds() if b.GetOtherAtomIdx(a.GetIdx()) in atoms)
    ]
    if doubles != 2 or len(sp3) != 1:
        return None
    return sp3[0]


def _stereoparent(bare):
    """('gonane' | 'morphinan', {atom of `bare`: locant}) when the bare ring system is that stereoparent, whose
    numbering is retained (P-26.2.2.2.1)."""
    if bare.GetNumAtoms() != 17 or bare.GetRingInfo().NumRings() < 4:
        return None
    if bare.HasSubstructMatch(_MORPHINAN):
        match = bare.GetSubstructMatch(_MORPHINAN)
        return "morphinan", {match[_MORPHINAN_INDEX[i]]: loc for i, loc in _MORPHINAN_LOCANTS.items()}
    if _is_saturated(bare) and all(a.GetAtomicNum() == 6 for a in bare.GetAtoms()):
        graph = adjacency(bare)
        found = _steroid_embeddings("gon", graph, set(range(bare.GetNumAtoms())))
        if len(found) == 1:
            return "gonane", {atom: position for position, atom in found[0].items()}
    return None


def _parent_rank(source, numbering):
    """Ring seniority of the parent hydride: the double bonds removed by 'hydro' prefixes still count (P-44.2)."""
    key = ring_seniority_key(source, range(source.GetNumAtoms()))
    return key[:-1] + (key[-1] - len(numbering.hydro_positions) // 2,)


def _route(bare):
    """'analogue': a saturated carbocycle, von Baeyer or spiro parent whose double and triple bonds are 'ene'/'yne'
    endings (P-31.1.6.1); 'mancude': a mancude parent, hydrogenation by 'hydro' prefixes (P-31.2.3.3.4)."""
    if any(a.GetIsAromatic() for a in bare.GetAtoms()):
        return "mancude"
    hetero = any(a.GetAtomicNum() != 6 for a in bare.GetAtoms())
    if bare.GetNumBonds() - bare.GetNumAtoms() + 1 == 1:
        return "analogue" if not hetero or bare.GetNumAtoms() > _LARGEST_HANTZSCH_WIDMAN_RING else "mancude"
    return "mancude" if is_hydro_fusion_system(bare, set(range(bare.GetNumAtoms()))) else "analogue"


def _aromatic_analogue(mol, atoms):
    source = Chem.RWMol()
    new_of = {}
    for old in sorted(atoms):
        atom = Chem.Atom(mol.GetAtomWithIdx(old).GetAtomicNum())
        atom.SetIsAromatic(True)
        new_of[old] = source.AddAtom(atom)
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms:
            source.AddBond(new_of[bond.GetBeginAtomIdx()], new_of[bond.GetEndAtomIdx()], Chem.BondType.AROMATIC)
    bare = source.GetMol()
    Chem.SanitizeMol(bare)
    return bare, {v: k for k, v in new_of.items()}


def build_amplificant(mol, atoms, free_atom=None, free_order=1):
    """The Amplificant of the ring system `atoms` of `mol`, or None when no parent hydride of P-26.2.2.2.1 describes it
    even with hydro prefixes, 'ene'/'yne' endings and skeletal replacement."""
    atoms = frozenset(atoms)
    added = None
    try:
        if free_order == 2 and free_atom in atoms:
            added = _added_hydrogen_ring(mol, atoms, free_atom)
            if added is None:
                return None
            bare, old_of = _aromatic_analogue(mol, atoms)
        else:
            bare, _, old_of = _bare_skeleton(mol, atoms)
    except Exception:
        return None
    stereoparent = _stereoparent(bare)
    if stereoparent is not None:
        name, locants = stereoparent
        return Amplificant(
            atoms=atoms,
            prefix=amplification_prefix(name),
            rank=ring_seniority_key(bare, range(bare.GetNumAtoms())) + (name,),
            mancude=any(a.GetIsAromatic() for a in bare.GetAtoms()),
            numberings=[Numbered({old_of[i]: AmpLoc(str(loc)) for i, loc in locants.items()})],
            added=added,
        )
    cyclomatic = bare.GetNumBonds() - bare.GetNumAtoms() + 1
    route = _route(bare)
    replaced, bonds, source = {}, [], bare
    if route == "analogue":
        heteroatoms = {i: a.GetSymbol() for i, a in enumerate(bare.GetAtoms()) if a.GetAtomicNum() != 6}
        if any(symbol not in _PREFIX for symbol in heteroatoms.values()):
            return None
        replaced = {old_of[i]: symbol for i, symbol in heteroatoms.items()}
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in atoms and b in atoms and bond.GetBondTypeAsDouble() > 1.0:
                bonds.append((a, b, int(bond.GetBondTypeAsDouble())))
        source = _saturated_analogue(bare, carbon=bool(heteroatoms))
    try:
        numberings = _numberings_of(source)
    except (UnsupportedStructure, ValueError, StopIteration, RuntimeError):
        return None
    if not numberings or any(n.added_positions for n in numberings):
        return None
    names = {_parent_name(n, source) for n in numberings}
    if len(names) != 1:
        return None
    name = names.pop()
    if ("spiro[" in name and cyclomatic > 2) or (route == "analogue" and any(n.unsat_key for n in numberings if getattr(n, "unsat_key", None))):
        return None
    best = min(n.pre_key for n in numberings)
    kept = [n for n in numberings if n.pre_key == best]
    amplificant = Amplificant(
        atoms=atoms,
        prefix=amplification_prefix(name),
        rank=_parent_rank(source, kept[0]) + (name,),
        mancude=route == "mancude"
        and not any(n.hydro_positions for n in kept)
        and (any(a.GetIsAromatic() for a in bare.GetAtoms()) or any(b.GetBondTypeAsDouble() == 2.0 for b in bare.GetBonds())),
        replaced=replaced,
        added=added,
        bonds=bonds,
        saturated_parent=route == "analogue",
    )
    for n in kept:
        local = {old_of[i]: AmpLoc(str(loc)) for i, loc in n.position_of.items()}
        ih = tuple(AmpLoc(str(p)) for p in getattr(n, "ih", ()))
        hydro = tuple(AmpLoc(str(p)) for p in n.hydro_positions)
        amplificant.numberings.append(Numbered(local, ih, hydro))
    return amplificant
