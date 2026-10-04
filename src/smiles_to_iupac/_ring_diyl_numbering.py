"""Numbering candidates and parent names for ring systems cited as multivalent groups (P-29.2, P-31.1.4):
monocycles (carbocycle, benzene, mancude/hydro/saturated heterocycles by retained and Hantzsch-Widman names,
P-22.2.2), ortho-fused arene chains, ortho-fused mancude heterocycles, von Baeyer and monospiro carbocycles.
"""

import re
from contextvars import ContextVar
from itertools import combinations

from rdkit import Chem

from . import _aromatic
from ._free_valence import valence_word
from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import (
    UnsupportedStructure,
    ring_bond_locants,
    von_baeyer_unsaturation_citations,
)
from ._fusion_numbering_hex import hex_numberings
from ._fusion_numbering_general import general_peripheral_numberings
from ._fusion_numbering_oriented import oriented_peripheral_numberings
from ._hetero_monocyclic import (
    _ROLE_SEQUENCES,
    saturated_five_membered_1_2_two_heteroatom_ring_name,
    saturated_five_membered_1_3_two_heteroatom_ring_name,
    saturated_ring_name,
    saturated_seven_membered_1_2_two_heteroatom_ring_name,
    saturated_seven_membered_1_3_two_heteroatom_ring_name,
    saturated_seven_membered_1_4_two_heteroatom_ring_name,
    saturated_two_heteroatom_1_4_ring_name,
)
from ._numerals import alkane_name, alkyl_name, multiplying_prefix
from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates
from ._spiro import find_monospiro_atom, iter_monospiro_numberings
from ._unsaturated import _unsaturation_suffix_from_citations
from ._common import multiplied_word

_SENIORITY = ["F", "Cl", "Br", "I", "O", "S", "Se", "Te", "N", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B", "Al", "Ga", "In", "Tl"]
_RANK = {e: i for i, e in enumerate(_SENIORITY)}
_PREFIX = {
    "O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha", "As": "arsa", "Sb": "stiba",
    "Bi": "bisma", "Si": "sila", "Ge": "germa", "Sn": "stanna", "Pb": "plumba", "B": "bora", "Al": "aluma",
    "Ga": "galla", "In": "indiga", "Tl": "thalla",
}
_CHALCOGENS = {"O", "S", "Se", "Te"}
_NO_DOUBLE_BOND = _CHALCOGENS | {"F", "Cl", "Br", "I"}
_SIX_A = {"O", "S", "Se", "Te", "Bi"}
_SIX_B = {"N", "Si", "Ge", "Sn", "Pb"}
_STEMS_UNSAT = {3: "irene", 4: "ete", 5: "ole", 7: "epine", 8: "ocine", 9: "onine", 10: "ecine"}
_STEMS_SAT = {3: "irane", 4: "etane", 5: "olane", 7: "epane", 8: "ocane", 9: "onane", 10: "ecane"}
_STEMS_SAT_N = {3: "iridine", 4: "etidine", 5: "olidine"}
_TWO_HETERO_SATURATED = {
    (6, 1, 4): saturated_two_heteroatom_1_4_ring_name,
    (5, 1, 3): saturated_five_membered_1_3_two_heteroatom_ring_name,
    (5, 1, 2): saturated_five_membered_1_2_two_heteroatom_ring_name,
    (7, 1, 4): saturated_seven_membered_1_4_two_heteroatom_ring_name,
    (7, 1, 3): saturated_seven_membered_1_3_two_heteroatom_ring_name,
    (7, 1, 2): saturated_seven_membered_1_2_two_heteroatom_ring_name,
}
_MANCUDE_RETAINED = {
    tuple(role[0] for role in roles): name[3:] if name.startswith("1H-") else name
    for name, roles in _ROLE_SEQUENCES.items()
    if name != "phosphinine"
}
_MANCUDE_RETAINED.update({(e, "C", "C", "C", "C", "C"): n for e, n in (("O", "pyran"), ("S", "thiopyran"), ("Se", "selenopyran"), ("Te", "telluropyran"))})
_MANCUDE_RETAINED[("N", "N", "N", "N", "C")] = "tetrazole"
_LACKS_REGULAR_NUMBERING = ("acridine", "carbazole", "xanthene", "purine", "anthracene", "phenanthrene")


class FusionLocant(float):
    """A lettered fusion-atom locant such as 4a: orders between 4 and 5, prints as written."""

    def __new__(cls, number, letters):
        obj = super().__new__(cls, number + sum((ord(c) - 96) / 100 ** (i + 1) for i, c in enumerate(letters)))
        obj.text = f"{number}{letters}"
        return obj

    __str__ = __repr__ = lambda self: self.text
    __format__ = lambda self, spec: self.text


def _locant(text):
    match = re.fullmatch(r"(\d+)([a-z]*)", text)
    if match is None:
        return None
    return FusionLocant(int(match.group(1)), match.group(2)) if match.group(2) else int(match.group(1))


SUFFIX_ATOMS = ContextVar("SUFFIX_ATOMS", default=frozenset())


class Numbering:
    def __init__(self, position_of, text, pre_key=(), unsat_key=(), ih=()):
        self.position_of = position_of
        self.text = text
        self.pre_key = pre_key
        self.unsat_key = unsat_key
        self.ih = ih
        self.ih_positions = ih
        self.hydro = ()
        self.added = ()
        self.fully_hydro = False
        self.stem = None


def _walks(ring_order):
    n = len(ring_order)
    for base in (ring_order, list(reversed(ring_order))):
        for start in range(n):
            yield base[start:] + base[:start]


def _yl(valence):
    return valence_word(valence)


def _tail(stem, locants, valence, suffix="yl"):
    if suffix == "carboxylate":
        return f"{stem}-{_locs(locants)}-carboxylate"
    if suffix == "":
        return stem
    if suffix != "yl":
        word = multiplied_word(valence, suffix)
        base = stem[:-1] if stem.endswith("e") and word[0] in "aeiouy" else stem
        return f"{base}-{_locs(locants)}-{word}"
    if valence == 1:
        base = stem[:-1] if stem.endswith("e") else stem
        return f"{base}-{locants[0]}-yl"
    return f"{stem}-{_locs(locants)}-{_yl(valence)}"


def _perfect_matching(atoms, adj):
    atoms = set(atoms)
    if not atoms:
        return True
    first = min(atoms)
    return any(
        _perfect_matching(atoms - {first, other}, adj) for other in adj[first] if other in atoms
    )


_RETAINED_NUMBERINGS = {
    "carbazole": (("1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "9a"), {"9": "N"}, (("4a", "9a"), ("4b", "8a"))),
    "fluorene": (("1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "9a"), {}, (("4a", "9a"), ("4b", "8a"))),
    "anthracene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {}, (("4a", "9a"), ("10a", "8a"))),
    "phenanthrene": (("1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "10", "10a"), {}, (("4a", "10a"), ("4b", "8a"))),
    "acridine": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "N"}, (("4a", "9a"), ("10a", "8a"))),
    "xanthene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "O"}, (("4a", "9a"), ("10a", "8a"))),
    "thioxanthene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "S"}, (("4a", "9a"), ("10a", "8a"))),
    "cyclopenta[a]phenanthrene": (
        ("1", "2", "3", "4", "5", "6", "7", "8", "14", "15", "16", "17", "13", "12", "11", "9", "10"),
        {},
        (("5", "10"), ("8", "9"), ("13", "14")),
    ),
    "purine": (("1", "2", "3", "4", "5", "6", "7", "8", "9"), {"1": "N", "3": "N", "7": "N", "9": "N"}, (("4", "5"),)),
    "fluoranthene": (
        ("1", "2", "3", "3a", "4", "5", "6", "6a", "6b", "7", "8", "9", "10", "10a", "10b", "3a1"),
        {},
        (("3a1", "3a"), ("3a1", "6a"), ("3a1", "10b"), ("6b", "10a")),
        15,
    ),
    "acenaphthylene": (
        ("1", "2", "2a", "3", "4", "5", "5a", "6", "7", "8", "8a", "2a1"),
        {},
        (("2a1", "2a"), ("2a1", "5a"), ("2a1", "8a")),
        11,
    ),
}


def _retained_numberings(bare, stem):
    entry = _RETAINED_NUMBERINGS[stem]
    locants, hetero, interior = entry[:3]
    index = {loc: i for i, loc in enumerate(locants)}
    ring_length = entry[3] if len(entry) == 4 else len(locants)
    edges = [(i, (i + 1) % ring_length) for i in range(ring_length)] + [(index[a], index[b]) for a, b in interior]
    if stem == "purine":
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (3, 8), (8, 7), (7, 6), (6, 4)]
    if bare.GetNumBonds() != len(edges):
        return []
    return [{atom: locants[i] for i, atom in enumerate(match)} for match in _template_matches(bare, locants, hetero, edges)]


def _template_matches(bare, locants, hetero, edges):
    symbol = {"N": 7, "O": 8, "S": 16}
    adjacency = {a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors()} for a in bare.GetAtoms()}
    wanted = [symbol.get(hetero.get(loc), 6) for loc in locants]
    neighbors = {i: set() for i in range(len(locants))}
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    if bare.GetNumAtoms() != len(locants):
        return []
    results = []

    def extend(assigned):
        i = len(assigned)
        if i == len(locants):
            results.append(tuple(assigned))
            return
        for atom in range(bare.GetNumAtoms()):
            if atom in assigned or bare.GetAtomWithIdx(atom).GetAtomicNum() != wanted[i]:
                continue
            if all((assigned[j] in adjacency[atom]) for j in neighbors[i] if j < i):
                extend(assigned + [atom])

    extend([])
    return results


def _split_hydrogen(position_of, adj, can_hold, saturated, oxo_all, oxo_suffix, ih_count):
    """(indicated, added, hydro) locant tuples for a ring whose `saturated` atoms include the ring atoms
    bearing oxo groups, following P-58.2.3: indicated hydrogen goes to the suffix atoms when it can
    accommodate them all, else the lowest valid position first; the suffix groups left over take added
    hydrogen (P-58.2.2). None when no split exists."""
    ordered = sorted(saturated, key=lambda a: position_of[a])
    if ih_count > len(ordered):
        return None
    loc = lambda atoms: tuple(sorted(position_of[a] for a in atoms))
    groups = len(oxo_suffix)

    def candidate_key(ih):
        if not groups:
            return (loc(ih),)
        valid = 0 if _perfect_matching(set(can_hold) - set(ih), adj) else 1
        covered = sum(a in oxo_suffix for a in ih)
        if ih_count >= groups:
            return (valid, -covered, loc(ih))
        return (valid, loc(ih), -covered)

    for ih in sorted(combinations(ordered, ih_count), key=candidate_key):
        free = [a for a in ordered if a not in ih and a not in oxo_suffix]
        base = set(can_hold) - set(ih) - set(oxo_suffix)
        for k in range(len(free) + 1):
            if (len(free) - k) % 2:
                continue
            found = next(
                (
                    added
                    for added in combinations(free, k)
                    if not oxo_suffix or _perfect_matching(base - set(added), adj)
                ),
                None,
            )
            if found is not None:
                return loc(ih), loc(found), loc(a for a in free if a not in found)
            if not oxo_suffix:
                break
    return None


def _tail_added(stem, locants, valence, suffix, added):
    if not added or suffix in ("", "carboxylate"):
        return _tail(stem, locants, valence, suffix)
    text = ",".join(f"{p}H" for p in added)
    if suffix == "yl":
        base = stem[:-1] if stem.endswith("e") and valence == 1 else stem
        return f"{base}-{_locs(locants)}({text})-{_yl(valence)}"
    word = multiplied_word(valence, suffix)
    base = stem[:-1] if stem.endswith("e") and word[0] in "aeiouy" else stem
    return f"{base}-{_locs(locants)}({text})-{word}"


def _ring_size(stem):
    from ._numerals import alkane_name as _an

    return next(n for n in range(3, 60) if _an(n) == stem[len("cyclo"):])


def _locs(locants):
    return ",".join(str(loc) for loc in sorted(locants))


def _hantzsch_widman(elements, saturated):
    size = len(elements)
    hetero = [e for e in elements if e != "C"]
    if size < 3 or size > 10 or not hetero or any(e not in _PREFIX for e in hetero):
        return None
    counts = {}
    for e in hetero:
        counts[e] = counts.get(e, 0) + 1
    ordered = sorted(counts, key=lambda e: _RANK[e])
    last = ordered[-1]
    if size == 6:
        group = "A" if last in _SIX_A else "B" if last in _SIX_B else "C"
        stem = ("ane" if group == "A" else "inane") if saturated else ("ine" if group != "C" else "inine")
    elif saturated:
        stem = _STEMS_SAT_N[size] if size in _STEMS_SAT_N and "N" in counts else _STEMS_SAT[size]
    else:
        stem = "irine" if size == 3 and set(counts) == {"N"} else _STEMS_UNSAT[size]
    pieces = []
    for e in ordered:
        mult = "" if counts[e] == 1 else multiplying_prefix(counts[e])
        pieces.append(mult + _PREFIX[e])
    text = ""
    for piece in pieces + [stem]:
        if text and piece[0] in "aeiou" and text.endswith("a"):
            text = text[:-1]
        text += piece
    return text


def _max_matching(atoms_can_hold, ring_order):
    n = len(ring_order)
    can = [a in atoms_can_hold for a in ring_order]
    edges = [(i, (i + 1) % n) for i in range(n) if can[i] and can[(i + 1) % n]]
    best = 0
    for r in range(len(edges), 0, -1):
        for combo in combinations(edges, r):
            used = [x for edge in combo for x in edge]
            if len(set(used)) == len(used):
                return r
    return best


def _hetero_monocycle(mol, ring_order, attached):
    atoms = [mol.GetAtomWithIdx(i) for i in ring_order]
    size = len(atoms)
    ring_set = set(ring_order)
    sym = {a.GetIdx(): a.GetSymbol() for a in atoms}
    if any(a.GetFormalCharge() != 0 for a in atoms):
        raise UnsupportedStructure("a charged ring atom is not supported as a diyl yet")
    aromatic = all(a.GetIsAromatic() for a in atoms)
    ring_double = {
        i
        for b in mol.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
        for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    }
    for a in atoms:
        for n in a.GetNeighbors():
            exocyclic_double = (
                n.GetIdx() not in ring_set and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
            )
            if exocyclic_double and (aromatic or ring_double) and a.GetAtomicNum() != 6:
                raise UnsupportedStructure("a ring atom with an exocyclic double bond is not supported as a diyl yet")
    if not aromatic and any(
        b.GetBondTypeAsDouble() not in (1.0, 2.0)
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
    ):
        raise UnsupportedStructure("a ring triple bond is not supported as a diyl yet")
    can_hold = {i for i in ring_order if sym[i] not in _NO_DOUBLE_BOND}
    oxo_all = _exocyclic_oxo(mol, ring_set)
    oxo_suffix = oxo_all & SUFFIX_ATOMS.get()
    if aromatic:
        in_double = {
            i
            for i in can_hold
            if i not in oxo_all
            and not (sym[i] != "C" and (mol.GetAtomWithIdx(i).GetTotalNumHs() > 0 or mol.GetAtomWithIdx(i).GetDegree() == 3))
        }
    else:
        in_double = ring_double
    saturated_atoms = can_hold - in_double
    oxo_suffix = oxo_suffix | {a for a in SUFFIX_ATOMS.get() & saturated_atoms if sym[a] != "C"}
    matching = _max_matching(can_hold, ring_order)
    mancude_sp3 = len(can_hold) - 2 * matching
    fully_saturated = not ring_double and not aromatic
    if not fully_saturated and not oxo_all:
        extra = len(saturated_atoms) - mancude_sp3
        if extra < 0 or extra % 2:
            raise UnsupportedStructure("this partly hydrogenated ring is not expressible with hydro prefixes yet")
    results = []
    best_pre = None
    senior = min(_RANK.get(sym[i], 99) for i in ring_order if sym[i] != "C")

    def starts_at_senior(walk):
        return _RANK.get(sym[walk[0]], 99) == senior

    for walk in _walks(ring_order):
        if not starts_at_senior(walk):
            continue
        hetero = [(i + 1, _RANK.get(sym[a], 99)) for i, a in enumerate(walk) if sym[a] != "C"]
        pre = (tuple(p for p, _ in hetero), tuple(r for _, r in hetero))
        if best_pre is None or pre < best_pre:
            best_pre = pre
    for walk in _walks(ring_order):
        if not starts_at_senior(walk):
            continue
        hetero = [(i + 1, _RANK.get(sym[a], 99)) for i, a in enumerate(walk) if sym[a] != "C"]
        pre = (tuple(p for p, _ in hetero), tuple(r for _, r in hetero))
        if pre != best_pre:
            continue
        position_of = {a: i + 1 for i, a in enumerate(walk)}
        elements = tuple(sym[a] for a in walk)
        if fully_saturated:
            stem = _saturated_name(elements, hetero)
            results.append((position_of, stem, (), (), ()))
            continue
        stem = _MANCUDE_RETAINED.get(elements)
        if stem is None:
            stem = _hantzsch_widman(elements, saturated=False)
            if stem is None:
                raise UnsupportedStructure("this heteromonocycle has no supported mancude parent name yet")
            stem = _with_hetero_locants(stem, elements, hetero)
        if oxo_all or oxo_suffix:
            adj = {a: {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in ring_set} for a in ring_order}
            split = _split_hydrogen(position_of, adj, can_hold, saturated_atoms, oxo_all, oxo_suffix, mancude_sp3)
            if split is None:
                continue
            ih, added, hydro = split
        else:
            sat_pos = sorted(position_of[a] for a in saturated_atoms)
            ih = tuple(sat_pos[:mancude_sp3])
            hydro = tuple(sat_pos[mancude_sp3:])
            added = ()
        results.append((position_of, stem, ih, hydro, added))
    return best_pre, results


def _saturated_name(elements, hetero):
    positions = [p for p, _ in hetero]
    size = len(elements)
    names = [e for e in elements if e != "C"]
    if len(names) == 1:
        stem = saturated_ring_name(names[0], size)
    elif len(names) == 2 and (size, positions[0], positions[1]) in _TWO_HETERO_SATURATED:
        stem = _TWO_HETERO_SATURATED[(size, positions[0], positions[1])](names)
    else:
        stem = None
    if stem is None:
        stem = _hantzsch_widman(elements, saturated=True)
    if stem is None:
        stem = _replacement_cycloalkane_name(elements)
    if stem is None:
        raise UnsupportedStructure("this saturated heteromonocycle has no supported name as a diyl yet")
    return _with_hetero_locants(stem, elements, hetero) if stem[0].isalpha() else stem


def _replacement_cycloalkane_name(elements):
    """Skeletal replacement name of a saturated ring with no Hantzsch-Widman name (P-22.2.3)."""
    by_element = {}
    for position, element in enumerate(elements, start=1):
        if element != "C":
            by_element.setdefault(element, []).append(position)
    if not by_element or any(e not in _PREFIX for e in by_element):
        return None
    pieces = []
    for element in sorted(by_element, key=lambda e: _RANK[e]):
        locants = by_element[element]
        mult = "" if len(locants) == 1 else multiplying_prefix(len(locants))
        pieces.append(f"{','.join(map(str, locants))}-{mult}{_PREFIX[element]}")
    return "-".join(pieces) + "cyclo" + alkane_name(len(elements))


_RETAINED_WITHOUT_LOCANTS = {"piperazine", "morpholine", "thiomorpholine", "imidazolidine", "pyrazolidine"}


def _with_hetero_locants(stem, elements, hetero):
    names = [e for e in elements if e != "C"]
    if stem in _RETAINED_WITHOUT_LOCANTS or len(names) < 2 or elements.count("C") <= 1:
        return stem
    pairs = sorted(zip(names, [p for p, _ in hetero]), key=lambda ep: (_RANK[ep[0]], ep[1]))
    return ",".join(str(p) for _, p in pairs) + "-" + stem


def monocycle_numberings(mol, ring_order, attached, valence, ene_bonds_getter=None):
    ring_set = set(ring_order)
    atoms = [mol.GetAtomWithIdx(i) for i in ring_order]
    size = len(ring_order)
    if all(a.GetIsAromatic() and a.GetAtomicNum() == 6 for a in atoms) and size == 6:
        out = []
        for walk in _walks(ring_order):
            position_of = {a: i + 1 for i, a in enumerate(walk)}
            out.append(Numbering(position_of, _benzene_text))
        return out
    if all(a.GetAtomicNum() == 6 for a in atoms):
        if any(a.GetIsAromatic() for a in atoms):
            raise UnsupportedStructure("an aromatic carbocycle other than benzene is not supported as a diyl yet")
        ring_bonds = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble())
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set and b.GetBondTypeAsDouble() != 1.0
        ]
        if any(order != 2.0 for _, _, order in ring_bonds):
            raise UnsupportedStructure("a ring triple bond is not supported as a diyl yet")
        stem = "cyclo" + alkane_name(size)
        out = []
        for walk in _walks(ring_order):
            position_of = {a: i + 1 for i, a in enumerate(walk)}
            ene = tuple(ring_bond_locants(position_of, ring_bonds, size)[0])
            out.append(Numbering(position_of, _carbocycle_text(stem, ene), unsat_key=ene))
        return out
    best_pre, results = _hetero_monocycle(mol, ring_order, attached)
    out = []
    for position_of, stem, ih, hydro, added in results:
        out.append(
            Numbering(
                position_of,
                _hetero_text(stem, ih, hydro, added),
                pre_key=best_pre + (ih,),
                unsat_key=(added, hydro),
                ih=ih,
            )
        )
    return out


def _benzene_text(locants, valence, substituted=frozenset(), suffix="yl"):
    loc = _locs(locants)
    if suffix == "carboxylate":
        return "benzoate"
    if valence == 1:
        return "phenyl"
    return f"{loc}-phenylene" if valence == 2 else f"benzene-{loc}-{_yl(valence)}"


def _carbocycle_text(stem, ene):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        if suffix == "carboxylate" and not ene and not substituted:
            return f"{stem}carboxylate"
        if not ene:
            return "cyclo" + alkyl_name(_ring_size(stem)) if valence == 1 and suffix == "yl" else _tail(stem, locants, valence, suffix)
        base = stem[:-3] + ("a" if len(ene) > 1 else "")
        return _tail(f"{base}-{','.join(map(str, ene))}-{multiplied_word(len(ene), 'ene')}", locants, valence, suffix)

    return text


def _hetero_text(stem, ih, hydro, added=()):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        ih_text = ",".join(f"{p}H" for p in ih) + "-" if ih else ""
        rest = ih_text + _tail_added(stem, locants, valence, suffix, added)
        hydro_text = f"{_locs(hydro)}-{multiplied_word(len(hydro), 'hydro')}" + ("-" if rest[0].isdigit() else "") if hydro else ""
        return hydro_text + rest

    return text



def _bare_skeleton(mol, skeleton_atoms, mancude=False):
    """Skeleton-only copy of the ring system; `mancude` makes every atom/bond aromatic (the mancude parent)."""
    order = sorted(skeleton_atoms)
    rw = Chem.RWMol()
    new_of = {}
    for old in order:
        atom = mol.GetAtomWithIdx(old)
        copy = Chem.Atom(atom.GetAtomicNum())
        copy.SetIsAromatic(True if mancude else atom.GetIsAromatic())
        copy.SetFormalCharge(atom.GetFormalCharge())
        new_of[old] = rw.AddAtom(copy)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in new_of and b in new_of:
            rw.AddBond(new_of[a], new_of[b], Chem.BondType.AROMATIC if mancude else bond.GetBondType())
    bare = rw.GetMol()
    for old in order:
        atom = mol.GetAtomWithIdx(old)
        if not mancude and atom.GetIsAromatic() and atom.GetAtomicNum() != 6:
            heavy_inside = sum(1 for n in atom.GetNeighbors() if n.GetIdx() in new_of)
            if atom.GetTotalNumHs() > 0 or atom.GetDegree() > heavy_inside or (atom.GetAtomicNum() == 7 and atom.GetDegree() == 3):
                bare.GetAtomWithIdx(new_of[old]).SetNumExplicitHs(1)
    if not mancude:
        Chem.SanitizeMol(bare)
    return bare, new_of, {v: k for k, v in new_of.items()}


ANION_SUFFIX = ContextVar("anion_suffix", default=frozenset())
_TRIVALENT_RING_HETERO = {7, 5, 13, 15, 31, 33, 49, 51, 81, 83}


def _mancude_candidates(mol, skeleton_atoms, sp3):
    """Skeleton-only mancude parents: every atom aromatic, with one explicit [nH] or one sp3 CH2 where needed."""
    bare, new_of, old_of = _bare_skeleton(mol, skeleton_atoms, mancude=True)
    nitrogens = [a for a in sorted(skeleton_atoms) if mol.GetAtomWithIdx(a).GetAtomicNum() in _TRIVALENT_RING_HETERO]
    attempts = [("none", None)] + [("nh", a) for a in nitrogens] + [
        ("ch2", a)
        for a in sorted(sp3, key=lambda a: (mol.GetRingInfo().NumAtomRings(a) > 1, a))
        if mol.GetAtomWithIdx(a).GetAtomicNum() == 6
    ]
    for mode, atom_idx in attempts:
        trial = Chem.RWMol(bare)
        if mode == "nh":
            trial.GetAtomWithIdx(new_of[atom_idx]).SetNumExplicitHs(1)
        elif mode == "ch2":
            target = new_of[atom_idx]
            trial.GetAtomWithIdx(target).SetIsAromatic(False)
            for bond in list(trial.GetAtomWithIdx(target).GetBonds()):
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
        try:
            Chem.SanitizeMol(trial)
        except Exception:
            continue
        yield trial.GetMol(), new_of, old_of


def _sanitized_mancude(mol, skeleton_atoms, sp3):
    for candidate in _mancude_candidates(mol, skeleton_atoms, sp3):
        return candidate
    raise UnsupportedStructure("the mancude parent of this ring system could not be built")


def _named_mancude(mol, skeleton_atoms, sp3):
    from .core import smiles_to_iupac

    for bare, new_of, old_of in _mancude_candidates(mol, skeleton_atoms, sp3):
        try:
            parent = smiles_to_iupac(Chem.MolToSmiles(bare))
        except UnsupportedStructure:
            continue
        if not re.search(r"cyclo\[|spiro\[|\d-hydro|\d-(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?ene$", parent):
            return bare, old_of, parent
    raise UnsupportedStructure("the mancude parent of this ring system has no fusion name")


def _hydro_text(hydro):
    return f"{_locs(hydro)}-{multiplied_word(len(hydro), 'hydro')}" if hydro else ""


def _vb_text(parent, ene_citations, hetero_prefix=""):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        stem = parent
        if ene_citations:
            body, needs_a = _unsaturation_suffix_from_citations(ene_citations, [])
            stem = parent[:-3] + ("a" if needs_a else "") + "-" + body
        return _tail(hetero_prefix + stem, locants, valence, suffix)

    return text


def _replacement_prefix(bare, order_new):
    """('2-oxa-5-aza' style text, hetero locants, element ranks) for a von Baeyer numbering."""
    by_element = {}
    for i, atom_idx in enumerate(order_new):
        symbol = bare.GetAtomWithIdx(atom_idx).GetSymbol()
        if symbol != "C":
            by_element.setdefault(symbol, []).append(i + 1)
    if not by_element:
        return "", (), ()
    if any(e not in _PREFIX for e in by_element):
        raise UnsupportedStructure("an unsupported skeletal heteroatom in a von Baeyer/spiro skeleton")
    pieces = []
    locants = []
    ranks = []
    for element in sorted(by_element, key=lambda e: _RANK[e]):
        locs = sorted(by_element[element])
        mult = "" if len(locs) == 1 else multiplying_prefix(len(locs))
        pieces.append(f"{','.join(map(str, locs))}-{mult}{_PREFIX[element]}")
        locants.extend(locs)
        ranks.extend([_RANK[element]] * len(locs))
    pairs = sorted(zip(locants, ranks))
    return "-".join(pieces), tuple(p for p, _ in pairs), tuple(r for _, r in pairs)


def _von_baeyer(mol, skeleton_atoms):
    bare, new_of, old_of = _bare_skeleton(mol, skeleton_atoms)
    if any(a.GetIsAromatic() for a in bare.GetAtoms()):
        raise UnsupportedStructure("an aromatic von Baeyer/spiro skeleton is not supported as a diyl yet")
    bonds_new = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble())
        for b in bare.GetBonds()
        if b.GetBondTypeAsDouble() != 1.0
    ]
    if any(order not in (2.0, 3.0) for _, _, order in bonds_new):
        raise UnsupportedStructure("an unsupported ring bond order in a polycyclic diyl")
    ring_count = bare.GetNumBonds() - bare.GetNumAtoms() + 1
    candidates = []
    spiro = find_monospiro_atom(bare) if ring_count == 2 else None
    if spiro is not None:
        for parent, order in iter_monospiro_numberings(bare, spiro):
            candidates.append((parent, order, ()))
    elif ring_count == 2:
        core = find_bicyclic_core(bare)
        if core is None:
            raise UnsupportedStructure("this bicyclic skeleton is not supported as a diyl yet")
        parent = bicyclic_parent_name(core)
        for order in iter_bicyclic_numberings(core):
            candidates.append((parent, order, ()))
    else:
        core = find_polycyclic_core(bare, ring_count)
        if core is None:
            raise UnsupportedStructure("this polycyclic skeleton is not supported as a diyl yet")
        for order, parent, outer in iter_polycyclic_candidates(core, ring_count):
            candidates.append((parent, order, outer))
    if not candidates:
        raise UnsupportedStructure("this ring skeleton has no supported von Baeyer/spiro numbering as a diyl yet")
    out = []
    for parent, order, outer in candidates:
        position_new = {a: i + 1 for i, a in enumerate(order)}
        position_of = {old_of[a]: p for a, p in position_new.items()}
        prefix, hetero_locs, hetero_ranks = _replacement_prefix(bare, order)
        pre = tuple(outer) + (hetero_locs, hetero_ranks)
        if bonds_new:
            ene, yne, compound_count, primary, full = von_baeyer_unsaturation_citations(position_new, bonds_new)
            if yne:
                raise UnsupportedStructure("a ring triple bond is not supported as a diyl yet")
            unsat = (compound_count, tuple(sorted(primary)), tuple(sorted(full)))
            out.append(Numbering(position_of, _vb_text(parent, ene, prefix), pre_key=pre, unsat_key=unsat))
        else:
            out.append(Numbering(position_of, _vb_text(parent, [], prefix), pre_key=pre))
    return out


def _sp3_ring_atoms(mol, skeleton_atoms):
    """Skeleton atoms that are neither aromatic nor in a skeleton double bond."""
    sp3 = set()
    for a in skeleton_atoms:
        atom = mol.GetAtomWithIdx(a)
        if atom.GetIsAromatic():
            continue
        if any(
            mol.GetBondBetweenAtoms(a, n.GetIdx()).GetBondTypeAsDouble() != 1.0 and n.GetIdx() in skeleton_atoms
            for n in atom.GetNeighbors()
        ):
            continue
        sp3.add(a)
    return sp3


def _arene_chain(mol, graph, rings, skeleton_atoms):
    atom_rings = [tuple(r) for r in rings]
    ring_atom_sets = [set(r) for r in atom_rings]
    fusion_bond_idxs = set()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if sum(1 for s in ring_atom_sets if a in s and b in s) >= 2:
            fusion_bond_idxs.add(bond.GetIdx())
    membership = _aromatic._atom_ring_membership(atom_rings)
    if any(len(r) >= 3 for r in membership.values()):
        raise UnsupportedStructure("a peri-fused arene is not supported as a diyl yet")
    adj, pairs = _aromatic._ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol)
    ring_order = _aromatic._ring_path_order(adj, len(atom_rings))
    shape = _aromatic._classify_shape(graph, atom_rings, ring_order, pairs)
    parent = _aromatic._retained_chain_name(len(atom_rings), shape)
    if parent is None:
        return None
    if parent == "anthracene":
        candidates = _aromatic._anthracene_candidates(graph, ring_atom_sets, pairs, ring_order)
    elif parent == "phenanthrene":
        candidates = _aromatic._phenanthrene_candidates(graph, ring_atom_sets, pairs, ring_order)
    else:
        candidates = _aromatic._straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order)
    oxo_all = _exocyclic_oxo(mol, skeleton_atoms)
    oxo_suffix = oxo_all & SUFFIX_ATOMS.get()
    sp3 = _sp3_ring_atoms(mol, skeleton_atoms) | oxo_all
    adj, can_hold = _ring_graph(mol, skeleton_atoms)
    out = []
    for c in candidates:
        position_of = dict(c)
        if any(a not in position_of for a in sp3):
            raise UnsupportedStructure("a saturated ring-fusion atom needs lettered hydro locants, not supported yet")
        split = _split_hydrogen(position_of, adj, can_hold, sp3, oxo_all, oxo_suffix, 0)
        if split is None:
            continue
        _, added, hydro = split

        def text(locants, valence, substituted=frozenset(), suffix="yl", hydro=hydro, added=added):
            return _hydro_text(hydro) + _tail_added(parent, locants, valence, suffix, added)

        out.append(Numbering(position_of, text, unsat_key=(added, hydro)))
    if not out:
        raise UnsupportedStructure("no numbering of this partly hydrogenated arene fits its hydro/added hydrogen")
    return out


def _exocyclic_oxo(mol, skeleton_atoms):
    return {
        a
        for a in skeleton_atoms
        if mol.GetAtomWithIdx(a).GetAtomicNum() == 6
        and any(
            n.GetIdx() not in skeleton_atoms and mol.GetBondBetweenAtoms(a, n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in mol.GetAtomWithIdx(a).GetNeighbors()
        )
    }


def _ring_graph(mol, skeleton_atoms):
    adj = {a: {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in skeleton_atoms} for a in skeleton_atoms}
    can_hold = {a for a in skeleton_atoms if mol.GetAtomWithIdx(a).GetSymbol() not in _NO_DOUBLE_BOND}
    return adj, can_hold


def _fused_mancude(mol, skeleton_atoms):
    oxo_all = _exocyclic_oxo(mol, skeleton_atoms)
    ring_info = mol.GetRingInfo()
    sp3 = {
        a for a in _sp3_ring_atoms(mol, skeleton_atoms) | oxo_all if mol.GetAtomWithIdx(a).GetSymbol() not in _NO_DOUBLE_BOND
    }
    sp3 |= {a for a in ANION_SUFFIX.get() if a in skeleton_atoms and mol.GetAtomWithIdx(a).GetAtomicNum() == 6}
    fusion_hetero = {a for a in skeleton_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and ring_info.NumAtomRings(a) > 1}
    suffix_atoms = SUFFIX_ATOMS.get() & set(skeleton_atoms)
    oxo_suffix = (oxo_all & suffix_atoms) | {a for a in suffix_atoms & sp3 if ring_info.NumAtomRings(a) > 1}
    bare, old_of, parent = _named_mancude(mol, skeleton_atoms, sp3)
    match = re.match(r"^(\d+H(?:,\d+H)*)-(.*)$", parent)
    ih_count = len(match.group(1).split(",")) if match else 0
    stem = match.group(2) if match else parent
    if stem in _RETAINED_NUMBERINGS:
        numberings = _retained_numberings(bare, stem)
    elif stem == "cyclopenta[a]phenanthrene" and not any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in skeleton_atoms):
        raise UnsupportedStructure("a saturated cyclopenta[a]phenanthrene skeleton is a steroid parent hydride (P-101), not a hydro fusion name")
    elif parent in _LACKS_REGULAR_NUMBERING:
        raise UnsupportedStructure("this fused parent has a retained non-peripheral numbering not supported yet")
    else:
        all_six = all(len(r) == 6 for r in bare.GetRingInfo().AtomRings())
        if all_six:
            numberings = hex_numberings(bare)
        else:
            try:
                numberings = oriented_peripheral_numberings(bare, ignore_indicated=True)
            except UnsupportedStructure:
                if bare.GetRingInfo().NumRings() > 3:
                    raise UnsupportedStructure("the numbering of this larger fused system with a peri-fused or other-sized ring is not verified")
                numberings = general_peripheral_numberings(bare, ignore_indicated=True)
    if not numberings:
        raise UnsupportedStructure("this fused skeleton has no supported peripheral numbering as a diyl yet")
    adj, can_hold = _ring_graph(mol, skeleton_atoms)
    can_hold -= fusion_hetero
    out = []
    for numbering in numberings:
        position_of = {old_of[new]: _locant(loc) for new, loc in numbering.items() if _locant(loc) is not None}
        if any(a not in position_of for a in sp3):
            continue
        saturated = {
            a
            for a in skeleton_atoms
            if a in position_of
            and (
                (a in sp3 and a not in fusion_hetero)
                or (
                    mol.GetAtomWithIdx(a).GetIsAromatic()
                    and mol.GetAtomWithIdx(a).GetAtomicNum() != 6
                    and (
                        mol.GetAtomWithIdx(a).GetTotalNumHs() > 0
                        or (mol.GetAtomWithIdx(a).GetDegree() == 3 and ring_info.NumAtomRings(a) < 2)
                    )
                )
            )
        }
        accommodated = oxo_suffix | {a for a in suffix_atoms & saturated if mol.GetAtomWithIdx(a).GetAtomicNum() != 6}
        if ANION_SUFFIX.get():
            centers = ANION_SUFFIX.get() & saturated
            accommodated = set() if saturated <= ANION_SUFFIX.get() else set(accommodated) | centers
        split = _split_hydrogen(position_of, adj, can_hold, saturated, oxo_all, accommodated, ih_count)
        if split is None:
            continue
        ih, added, hydro = split
        ih_text = (",".join(f"{p}H" for p in ih) + "-") if ih else ""

        fully_saturated = bool(hydro) and can_hold <= saturated | set(accommodated)

        def text(locants, valence, substituted=frozenset(), suffix="yl", hydro=hydro, ih_text=ih_text, added=added, full=fully_saturated):
            rest = ih_text + _tail_added(stem, locants, valence, suffix, added)
            hydro_text = multiplied_word(len(hydro), "hydro") if full else _hydro_text(hydro)
            return hydro_text + ("-" if hydro and rest[0].isdigit() else "") + rest

        numbering = Numbering(position_of, text, pre_key=(ih,), unsat_key=(added, hydro), ih=ih)
        numbering.hydro, numbering.added, numbering.fully_hydro, numbering.stem = tuple(hydro), tuple(added), fully_saturated, stem
        out.append(numbering)
    if not out:
        raise UnsupportedStructure("no numbering of this partly hydrogenated fused system fits its hydro/indicated hydrogen")
    return out


def is_hydro_fusion_system(mol, skeleton_atoms):
    """A non-aromatic ortho-fused ring system with at least two rings of five or more members
    is named as a hydro derivative of its mancude fusion parent (P-31.1.4.2.4), not by von Baeyer."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= set(skeleton_atoms)]
    if len(rings) < 2 or sum(len(r) >= 5 for r in rings) < 2:
        return False
    if any(a.GetIsAromatic() for a in map(mol.GetAtomWithIdx, skeleton_atoms)):
        return False
    pairs = [(i, j) for i in range(len(rings)) for j in range(i + 1, len(rings)) if rings[i] & rings[j]]
    if len(pairs) != len(rings) - 1 or any(len(rings[i] & rings[j]) != 2 for i, j in pairs):
        return False
    return not any(sum(a in r for r in rings) > 2 for a in skeleton_atoms)


def system_numberings(mol, graph, rings, skeleton_atoms):
    if all(
        len(r) == 6 and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in r) for r in rings
    ) and len(rings) <= 3 and any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in skeleton_atoms):
        try:
            arene = _arene_chain(mol, graph, rings, skeleton_atoms)
        except UnsupportedStructure:
            arene = None
        if arene is not None:
            return arene
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in skeleton_atoms):
        return _fused_mancude(mol, skeleton_atoms)
    if is_hydro_fusion_system(mol, skeleton_atoms):
        try:
            return _fused_mancude(mol, skeleton_atoms)
        except UnsupportedStructure:
            pass
    return _von_baeyer(mol, skeleton_atoms)


def chain_numberings(mol, paths, valence):
    from ._common import name_from_substituents

    out = []
    for path in paths:
        path_set = set(path)
        for atom in path:
            for n in mol.GetAtomWithIdx(atom).GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom, n.GetIdx())
                if n.GetIdx() not in path_set and bond.GetBondTypeAsDouble() != 1.0 and n.GetAtomicNum() == 6:
                    raise UnsupportedStructure("an ylidene-type multiple bond off the polyol chain is not supported yet")
        for direction in (path, list(reversed(path))):
            position_of = {a: i + 1 for i, a in enumerate(direction)}
            ene, yne = [], []
            for i in range(len(direction) - 1):
                bond = mol.GetBondBetweenAtoms(direction[i], direction[i + 1])
                order = bond.GetBondTypeAsDouble()
                if order == 2.0:
                    ene.append(i + 1)
                elif order == 3.0:
                    yne.append(i + 1)
                elif order != 1.0:
                    raise UnsupportedStructure("an aromatic or unusual bond on the polyol chain")
            length = len(direction)

            def text(locants, valence, substituted=frozenset(), suffix="yl", ene=tuple(ene), yne=tuple(yne), length=length):
                if length == 1 and valence == 2:
                    return "methylene"
                if valence == 1 and not ene and not yne:
                    return alkyl_name(length) if locants[0] == 1 else name_from_substituents(length, [], [], "yl", own_locants=list(locants))
                return name_from_substituents(length, list(ene), list(yne), _yl(valence), own_locants=sorted(locants))

            out.append(Numbering(position_of, text, unsat_key=(tuple(ene), tuple(yne))))
    return out
