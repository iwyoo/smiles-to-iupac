"""Numbering candidates and parent names for ring systems cited as multivalent groups (P-29.2, P-31.1.4):
monocycles (carbocycle, benzene, mancude/hydro/saturated heterocycles by retained and Hantzsch-Widman names,
P-22.2.2), ortho-fused arene chains, ortho-fused mancude heterocycles, von Baeyer and monospiro carbocycles.
"""

import re
import networkx as nx
from contextvars import ContextVar
from itertools import combinations

from rdkit import Chem

from ._free_valence import valence_word
from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import (
    UnsupportedStructure,
    nonstandard_bonding,
    ring_bond_locants,
    von_baeyer_unsaturation_citations,
)
from ._fusion_name import FUSION_NAME_REQUIRED
from ._fusion_components import EXCEPTIONS, system_numberings as fusion_system_numberings
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
from ._polyspiro import (
    find_branched_polyspiro_hub,
    find_linear_polyspiro_chain,
    iter_branched_polyspiro_numberings,
    iter_linear_polyspiro_numberings,
)
from ._spiro import find_monospiro_atom, iter_monospiro_numberings
from ._unsaturated import _unsaturation_suffix_from_citations
from ._common import multiplied_word, reparsable_smiles, sanitize_probe

_SENIORITY = ["F", "Cl", "Br", "I", "O", "S", "Se", "Te", "N", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B", "Al", "Ga", "In", "Tl"]
_RANK = {e: i for i, e in enumerate(_SENIORITY)}
_PREFIX = {
    "F": "fluora", "Cl": "chlora", "Br": "broma", "I": "ioda",
    "O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha", "As": "arsa", "Sb": "stiba",
    "Bi": "bisma", "Si": "sila", "Ge": "germa", "Sn": "stanna", "Pb": "plumba", "B": "bora", "Al": "aluma",
    "Ga": "galla", "In": "indiga", "Tl": "thalla",
}
def _lambda_by_position(mol, position_of):
    return {p: n for a, p in position_of.items() if (n := nonstandard_bonding(mol.GetAtomWithIdx(a)))}


def _lambda_key(lam):
    """P-14.4(h), P-22.2.7.2: low locants to the atoms with nonstandard bonding numbers, the higher number first on a tie."""
    return tuple(sorted(lam)), tuple(-lam[p] for p in sorted(lam))


def _cite(locant, lam):
    return f"{locant}\u03bb{lam[locant]}" if locant in lam else str(locant)
_CHALCOGENS = {"O", "S", "Se", "Te"}
_NO_DOUBLE_BOND = _CHALCOGENS | {"F", "Cl", "Br", "I"}
_TETRAVALENT = {"C", "Si", "Ge", "Sn", "Pb"}
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


class FusionLocant(float):
    """A lettered fusion-atom locant such as 4a: orders between 4 and 5, prints as written."""

    def __new__(cls, number, letters, interior=""):
        value = number + sum((ord(c) - 96) / 100 ** (i + 1) for i, c in enumerate(letters))
        obj = super().__new__(cls, value + (int(interior) / 100 ** (len(letters) + 1) if interior else 0))
        obj.text = f"{number}{letters}{interior}"
        return obj

    __str__ = __repr__ = lambda self: self.text
    __format__ = lambda self, spec: self.text


def _locant(text):
    match = re.fullmatch(r"(\d+)([a-z]*)(\d*)", text)
    if match is None:
        return None
    return FusionLocant(int(match.group(1)), match.group(2), match.group(3)) if match.group(2) else int(match.group(1))


SUFFIX_ATOMS = ContextVar("SUFFIX_ATOMS", default=frozenset())


class Numbering:
    def __init__(self, position_of, text, pre_key=(), unsat_key=(), ih=()):
        self.position_of = position_of
        self.text = text
        self.pre_key = pre_key
        self.unsat_key = unsat_key
        self.lam_key = ()
        self.ih = ih
        self.ih_positions = ih
        self.hydro = ()
        self.added = ()
        self.fully_hydro = False
        self.stem = None
        self.parent_stem = None
        self.hydro_positions = ()
        self.added_positions = ()


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

    def added_needed(ih):
        free = [a for a in ordered if a not in ih and a not in oxo_suffix]
        base = set(can_hold) - set(ih) - set(oxo_suffix)
        for k in range(len(free) + 1):
            if (len(free) - k) % 2:
                continue
            if any(_perfect_matching(base - set(added), adj) for added in combinations(free, k)):
                return k
        return len(free) + 1

    def candidate_key(ih):
        valid = 0 if _perfect_matching(set(can_hold) - set(ih), adj) else 1
        if not groups:
            return (valid, loc(ih))
        # P-58.2.3.1.4: indicated hydrogen that makes added hydrogen unnecessary comes before a placement that needs it
        needed = added_needed(ih)
        covered = sum(a in oxo_suffix for a in ih)
        if ih_count >= groups:
            return (valid, needed, -covered, loc(ih))
        return (valid, needed, loc(ih), -covered)

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


def _perfect_matchings(vertices, adj, limit=2000):
    found = []

    def extend(left, chosen):
        if len(found) > limit:
            return
        if not left:
            found.append(frozenset(chosen))
            return
        first = min(left)
        for other in sorted(adj[first] & left):
            extend(left - {first, other}, chosen + [frozenset((first, other))])

    extend(frozenset(vertices), [])
    return found


def _delta_citation(mol, skeleton_atoms, position_of):
    """P-25.7.1.2: the localized double bonds (as the Greek capital delta locants) that tell this isomer from the other
    arrangements of the same ring atoms and substituents; () when the structure is the only one or delocalized."""
    ring = set(skeleton_atoms)
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return ()
    double = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in mol.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring
    }
    live = set().union(*double) if double else set()
    if len(live) < 4 or len(live) > 24:
        return ()
    adj = {a: {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in live} for a in live}
    matchings = _perfect_matchings(live, adj)
    if len(matchings) < 2 or len(matchings) > 2000:
        return ()
    ring_bonds = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds() if b.GetBeginAtomIdx() in live and b.GetEndAtomIdx() in live}

    def identity(matching):
        trial = Chem.RWMol(mol)
        for bond in ring_bonds:
            a, b = tuple(bond)
            trial.GetBondBetweenAtoms(a, b).SetBondType(Chem.BondType.DOUBLE if bond in matching else Chem.BondType.SINGLE)
        out = trial.GetMol()
        Chem.SanitizeMol(out, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
        return Chem.MolToSmiles(out)

    actual = frozenset(double & ring_bonds)
    own = identity(actual)
    others = [m for m in matchings if identity(m) != own]
    if not others:
        return ()
    order = sorted(live, key=lambda a: position_of[a])
    following = {a: order[i + 1] for i, a in enumerate(order[:-1])}

    def token(bond):
        a, b = sorted(bond, key=lambda x: position_of[x])
        return (position_of[a], str(position_of[a]) if following.get(a) == b else f"{position_of[a]}({position_of[b]})")

    bonds = sorted(actual, key=lambda bond: sorted(position_of[x] for x in bond))
    for size in range(1, len(bonds) + 1):
        for subset in combinations(bonds, size):
            if not any(set(subset) <= m for m in others):
                return tuple(sorted(token(bond) for bond in subset))
    return ()


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
        if mult.endswith("a") and _PREFIX[e][0] in "aeiou":
            mult = mult[:-1]
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
    triple_atoms = {
        i
        for b in mol.GetBonds()
        if b.GetBondTypeAsDouble() == 3.0 and b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
        for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    }
    ring_double |= triple_atoms
    heterones = _chalcogen_heterones(mol, ring_set) & SUFFIX_ATOMS.get() if not aromatic else set()
    for a in atoms:
        for n in a.GetNeighbors():
            exocyclic_double = (
                n.GetIdx() not in ring_set and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
            )
            if exocyclic_double and (aromatic or ring_double) and a.GetAtomicNum() != 6 and a.GetIdx() not in heterones:
                raise UnsupportedStructure("a ring atom with an exocyclic double bond is not supported as a diyl yet")
    if not aromatic and any(
        b.GetBondTypeAsDouble() not in (1.0, 2.0, 3.0)
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
    ):
        raise UnsupportedStructure("this ring bond is not supported as a diyl yet")
    can_hold = {i for i in ring_order if sym[i] not in _NO_DOUBLE_BOND or (sym[i] in ("S", "Se", "Te") and i in ring_double)} | heterones
    oxo_all = _exocyclic_oxo(mol, ring_set) | (heterones - ring_double)
    oxo_suffix = oxo_all & SUFFIX_ATOMS.get()
    if aromatic:
        in_double = {
            i
            for i in can_hold
            if i not in oxo_all
            and not (
                sym[i] not in _TETRAVALENT
                and (mol.GetAtomWithIdx(i).GetTotalNumHs() > 0 or mol.GetAtomWithIdx(i).GetDegree() == 3)
                and not mol.GetAtomWithIdx(i).HasProp("_ring_cation_centre")
            )
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
        position_of = {a: i + 1 for i, a in enumerate(walk)}
        lam = _lambda_by_position(mol, position_of)
        lam.update({position_of[a]: int(mol.GetAtomWithIdx(a).GetTotalValence()) for a in heterones})
        pre = (tuple(p for p, _ in hetero), tuple(r for _, r in hetero))
        if pre != best_pre:
            continue
        elements = tuple(sym[a] for a in walk)
        if fully_saturated:
            stem = _saturated_name(elements, hetero, lam)
            results.append((position_of, stem, (), (), (), ()))
            continue
        stem = _MANCUDE_RETAINED.get(elements)
        if stem is None:
            stem = _hantzsch_widman(elements, saturated=False)
            if stem is None:
                stem = _replacement_ene_name(mol, elements, walk, hetero, lam)
                if stem is None:
                    raise UnsupportedStructure("this heteromonocycle has no supported mancude parent name yet")
                results.append((position_of, stem, (), (), (), tuple(sorted(position_of[a] for a in triple_atoms))))
                continue
            stem = _with_hetero_locants(stem, elements, hetero, lam)
        elif lam:
            stem = _with_hetero_locants(re.sub(r"^\d+(?:,\d+)*-", "", stem), elements, hetero, lam)
        if oxo_all or oxo_suffix:
            adj = {a: {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in ring_set} for a in ring_order}
            split = _split_hydrogen(position_of, adj, can_hold, saturated_atoms, oxo_all, oxo_suffix, mancude_sp3)
            if split is None:
                continue
            ih, added, hydro = split
        else:
            sat_pos = sorted(position_of[a] for a in saturated_atoms)
            lambda_h = [
                position_of[a]
                for a in walk
                if position_of[a] in lam
                and mol.GetAtomWithIdx(a).GetTotalNumHs()
                and not any(b.GetBondTypeAsDouble() == 2.0 for b in mol.GetAtomWithIdx(a).GetBonds())
            ]
            ih = tuple(sorted(sat_pos[:mancude_sp3] + lambda_h))
            hydro = tuple(sat_pos[mancude_sp3:])
            added = ()
        results.append((position_of, stem, ih, hydro, added, tuple(sorted(position_of[a] for a in triple_atoms))))
    return best_pre, results


def _saturated_name(elements, hetero, lam=None):
    positions = [p for p, _ in hetero]
    size = len(elements)
    names = [e for e in elements if e != "C"]
    if len(names) == 1:
        stem = saturated_ring_name(names[0], size)
    elif not lam and len(names) == 2 and (size, positions[0], positions[1]) in _TWO_HETERO_SATURATED:
        stem = _TWO_HETERO_SATURATED[(size, positions[0], positions[1])](names)
    else:
        stem = None
    if stem is None:
        stem = _hantzsch_widman(elements, saturated=True)
    if stem is None:
        stem = _replacement_cycloalkane_name(elements, lam)
    if stem is None:
        raise UnsupportedStructure("this saturated heteromonocycle has no supported name as a diyl yet")
    return _with_hetero_locants(stem, elements, hetero, lam) if stem[0].isalpha() else stem


def _replacement_cycloalkane_name(elements, lam=None):
    """Skeletal replacement name of a saturated ring with no Hantzsch-Widman name (P-22.2.3)."""
    by_element = {}
    for position, element in enumerate(elements, start=1):
        if element != "C":
            by_element.setdefault(element, []).append(position)
    if not by_element or any(e not in _PREFIX for e in by_element):
        return None
    lam = lam or {}
    if not lam and len(by_element) == 1 and len(next(iter(by_element.values()))) == len(elements):
        (element,) = by_element
        return multiplying_prefix(len(elements)) + _PREFIX[element] + "cyclo" + alkane_name(len(elements))
    pieces = []
    for element in sorted(by_element, key=lambda e: _RANK[e]):
        locants = by_element[element]
        mult = "" if len(locants) == 1 else multiplying_prefix(len(locants))
        pieces.append(f"{','.join(_cite(p, lam) for p in locants)}-{mult}{_PREFIX[element]}")
    return "-".join(pieces) + "cyclo" + alkane_name(len(elements))


def _replacement_ene_name(mol, elements, walk, hetero, lam):
    """Replacement name of a ring too large for a Hantzsch-Widman name, its double bonds cited as 'ene' (P-22.2.3, P-31.1.4.2.4)."""
    if len(elements) < 11:
        return None
    base = _replacement_cycloalkane_name(elements, lam)
    if base is None:
        return None
    position_of = {a: i + 1 for i, a in enumerate(walk)}
    enes = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bond.GetBondTypeAsDouble() == 2.0 and a in position_of and b in position_of:
            low, high = sorted((position_of[a], position_of[b]))
            enes.append((low, str(low)) if high - low == 1 else (low, f"{low}({high})"))
    body, needs_a = _unsaturation_suffix_from_citations(enes, [])
    return base[:-3] + ("a" if needs_a else "") + "-" + body


_RETAINED_WITHOUT_LOCANTS = {"piperazine", "morpholine", "thiomorpholine", "imidazolidine", "pyrazolidine"}


def _with_hetero_locants(stem, elements, hetero, lam=None):
    lam = lam or {}
    names = [e for e in elements if e != "C"]
    if not lam and (
        stem in _RETAINED_WITHOUT_LOCANTS
        or len(names) < 2
        or (elements.count("C") <= 1 and len(set(names)) < 2)
        or (len(elements) == 3 and len(names) == 2)
    ):
        return stem
    pairs = sorted(zip(names, [p for p, _ in hetero]), key=lambda ep: (_RANK[ep[0]], ep[1]))
    return ",".join(_cite(p, lam) for _, p in pairs) + "-" + stem


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
    for position_of, stem, ih, hydro, added, dehydro in results:
        numbering = Numbering(
            position_of,
            _hetero_text(stem, ih, hydro, added, dehydro),
            pre_key=best_pre + (ih,),
            unsat_key=(added, tuple(sorted(hydro + dehydro))),
            ih=ih,
        )
        numbering.parent_stem, numbering.hydro_positions, numbering.added_positions = stem, tuple(hydro), tuple(added)
        numbering.lam_key = _lambda_key(_lambda_by_position(mol, position_of))
        out.append(numbering)
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


def _hetero_text(stem, ih, hydro, added=(), dehydro=()):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        # P-14.3.4.2(c): the position of the one heteroatom of an unsubstituted ring is not cited (thiacyclododecane)
        bare = stem
        if not substituted and not valence and suffix in ("", "yl") and not hydro and not added:
            bare = re.sub(r"^1-(?=[a-z]+acyclo)", "", stem)
        ih_text = ",".join(f"{p}H" for p in ih) + "-" if ih else ""
        rest = ih_text + _tail_added(bare, locants, valence, suffix, added)
        hydro_text = f"{_locs(hydro)}-{multiplied_word(len(hydro), 'hydro')}" + ("-" if rest[0].isdigit() else "") if hydro else ""
        if dehydro:
            hydro_text = f"{_locs(dehydro)}-{multiplied_word(len(dehydro), 'dehydro')}" + (
                "-" if hydro_text or rest[0].isdigit() else ""
            ) + hydro_text
        return hydro_text + rest

    return text



def _bare_skeleton(mol, skeleton_atoms, mancude=False):
    """Skeleton-only copy of the ring system; `mancude` makes every atom/bond aromatic (the mancude parent)."""
    order = sorted(skeleton_atoms)
    rw = Chem.RWMol()
    new_of = {}
    ring_halogens = {old for old in order if mancude and mol.GetAtomWithIdx(old).GetAtomicNum() in (9, 17, 35, 53)}
    for old in order:
        atom = mol.GetAtomWithIdx(old)
        copy = Chem.Atom(atom.GetAtomicNum())
        copy.SetIsAromatic(True if mancude and old not in ring_halogens else atom.GetIsAromatic())
        copy.SetFormalCharge(atom.GetFormalCharge())
        if old in ring_halogens:
            copy.SetNoImplicit(True)
            copy.SetNumExplicitHs(1)
        new_of[old] = rw.AddAtom(copy)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in new_of and b in new_of:
            single = a in ring_halogens or b in ring_halogens
            rw.AddBond(
                new_of[a], new_of[b], Chem.BondType.SINGLE if single else Chem.BondType.AROMATIC if mancude else bond.GetBondType()
            )
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
    carbons = [
        a
        for a in sorted(sp3, key=lambda a: (mol.GetRingInfo().NumAtomRings(a) > 1, a))
        if mol.GetAtomWithIdx(a).GetAtomicNum() in (6, 14, 32, 50, 82)
    ]
    attempts = (
        [((), ())]
        + [((a,), ()) for a in nitrogens]
        + [((), (a,)) for a in carbons]
        + [((), pair) for pair in combinations(carbons, 2)]
        + [(pair, ()) for pair in combinations(nitrogens, 2)]
        + [((n,), (c,)) for n in nitrogens for c in carbons]
    )
    holders = {
        a
        for a in skeleton_atoms
        if mol.GetAtomWithIdx(a).GetSymbol() not in _NO_DOUBLE_BOND
        and not (
            mol.GetAtomWithIdx(a).GetAtomicNum() in _TRIVALENT_RING_HETERO
            and mol.GetAtomWithIdx(a).GetDegree() == 3
            and not mol.GetAtomWithIdx(a).GetFormalCharge()
            and mol.GetRingInfo().NumAtomRings(a) > 1
        )
    }
    neighbours = {
        a: {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in holders} for a in holders
    }
    for hydrogenated, saturated in attempts:
        chosen = set(hydrogenated) | set(saturated)
        if not _perfect_matching(holders - chosen, neighbours):
            continue
        trial = Chem.RWMol(bare)
        for member in hydrogenated:
            trial.GetAtomWithIdx(new_of[member]).SetNumExplicitHs(1)
        for member in saturated:
            target = new_of[member]
            trial.GetAtomWithIdx(target).SetIsAromatic(False)
            for bond in list(trial.GetAtomWithIdx(target).GetBonds()):
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
            if trial.GetAtomWithIdx(target).GetAtomicNum() in (14, 32, 50, 82):
                trial.GetAtomWithIdx(target).SetNoImplicit(True)
                trial.GetAtomWithIdx(target).SetNumExplicitHs(4 - trial.GetAtomWithIdx(target).GetDegree())
        try:
            sanitize_probe(trial)
        except Exception:
            continue
        yield trial.GetMol(), new_of, old_of


_MANCUDE_IN_PROGRESS = ContextVar("mancude_in_progress", default=frozenset())


def _named_mancude(mol, skeleton_atoms, sp3):
    from .core import smiles_to_iupac

    for bare, new_of, old_of in _mancude_candidates(mol, skeleton_atoms, sp3):
        smiles = reparsable_smiles(bare)
        if smiles in _MANCUDE_IN_PROGRESS.get():
            continue
        token = _MANCUDE_IN_PROGRESS.set(_MANCUDE_IN_PROGRESS.get() | {smiles})
        fusion_token = FUSION_NAME_REQUIRED.set(True)
        try:
            parent = smiles_to_iupac(smiles)
        except UnsupportedStructure:
            continue
        finally:
            FUSION_NAME_REQUIRED.reset(fusion_token)
            _MANCUDE_IN_PROGRESS.reset(token)
        if not re.search(r"cyclo\[|spiro\[|\d-hydro|\d-(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?ene$", parent):
            return bare, old_of, parent
    raise UnsupportedStructure("the mancude parent of this ring system has no fusion name")


def _hydro_text(hydro):
    return f"{_locs(hydro)}-{multiplied_word(len(hydro), 'hydro')}" if hydro else ""


def _vb_text(parent, ene_citations, hetero_prefix="", yne_citations=()):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        stem = parent
        if ene_citations or yne_citations:
            body, needs_a = _unsaturation_suffix_from_citations(ene_citations, yne_citations)
            stem = parent[:-3] + ("a" if needs_a else "") + "-" + body
        return _tail(hetero_prefix + stem, locants, valence, suffix)

    return text


def _replacement_prefix(bare, order_new, lam_new=None):
    """('2-oxa-5-aza' style text, hetero locants, element ranks, lambda key) for a von Baeyer numbering; `lam_new` maps
    skeleton atoms to a nonstandard bonding number (P-15.4.1.3)."""
    lam = {i + 1: lam_new[atom_idx] for i, atom_idx in enumerate(order_new) if atom_idx in (lam_new or {})}
    by_element = {}
    for i, atom_idx in enumerate(order_new):
        symbol = bare.GetAtomWithIdx(atom_idx).GetSymbol()
        if symbol != "C":
            by_element.setdefault(symbol, []).append(i + 1)
    if not by_element:
        return "", (), (), ()
    if any(e not in _PREFIX for e in by_element):
        raise UnsupportedStructure("an unsupported skeletal heteroatom in a von Baeyer/spiro skeleton")
    pieces = []
    locants = []
    ranks = []
    for element in sorted(by_element, key=lambda e: _RANK[e]):
        locs = sorted(by_element[element])
        mult = "" if len(locs) == 1 else multiplying_prefix(len(locs))
        pieces.append(f"{','.join(_cite(p, lam) for p in locs)}-{mult}{_PREFIX[element]}")
        locants.extend(locs)
        ranks.extend([_RANK[element]] * len(locs))
    pairs = sorted(zip(locants, ranks))
    return "-".join(pieces), tuple(p for p, _ in pairs), tuple(r for _, r in pairs), _lambda_key(lam)


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
    elif (chain := find_linear_polyspiro_chain(bare)) is not None:
        for parent, order, spiros, _ in iter_linear_polyspiro_numberings(bare, chain):
            candidates.append((parent, order, ()))
    elif (hub := find_branched_polyspiro_hub(bare)) is not None:
        for parent, order, spiros, _ in iter_branched_polyspiro_numberings(bare, hub):
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
        lam_new = {a: n for a, old_idx in old_of.items() if (n := nonstandard_bonding(mol.GetAtomWithIdx(old_idx)))}
        prefix, hetero_locs, hetero_ranks, lam_key = _replacement_prefix(bare, order, lam_new)
        pre = tuple(outer) + (hetero_locs, hetero_ranks, lam_key)
        if bonds_new:
            ene, yne, compound_count, primary, full = von_baeyer_unsaturation_citations(position_new, bonds_new)
            unsat = (compound_count, tuple(sorted(primary)), tuple(sorted(full)), tuple(sorted(p for p, _ in ene)))
            out.append(Numbering(position_of, _vb_text(parent, ene, prefix, yne), pre_key=pre, unsat_key=unsat))
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


def _ring_halogens(mol, skeleton_atoms):
    """Halogen atoms bonded to two ring atoms: λ3 (or λ5) centres that take part in a ring double bond of the mancude
    parent (P-22.2.7.1)."""
    return {
        a
        for a in skeleton_atoms
        if mol.GetAtomWithIdx(a).GetSymbol() in ("Cl", "Br", "I")
        and sum(1 for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in skeleton_atoms) == 2
    }


def _ring_graph(mol, skeleton_atoms):
    adj = {a: {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in skeleton_atoms} for a in skeleton_atoms}
    can_hold = {a for a in skeleton_atoms if mol.GetAtomWithIdx(a).GetSymbol() not in _NO_DOUBLE_BOND}
    return adj, can_hold | _ring_halogens(mol, skeleton_atoms)


def _replacement_numberings(bare):
    """Numberings of a system named by skeletal replacement on a hydrocarbon fusion name (P-25.5.1): the hydrocarbon's
    numbering is kept and the heteroatoms take the lowest locants it permits."""
    from ._fusion_name import _ReplacementRoot, fusion_name

    try:
        _, root = fusion_name(bare)
    except UnsupportedStructure:
        return None
    return list(root.numberings) if isinstance(root, _ReplacementRoot) else None


def _cite_lambda(stem, bonding):
    """(stem, text) with the λ marks joined to the heteroatom locants a fusion name cites in front ('1λ6,2-benzothiazole') and
    the others cited loose before the name (P-25.6)."""
    from ._fused_numbering import _locant_key

    leading = re.match(r"(\d+[a-z]?(?:,\d+[a-z]?)*)-(.+)", stem)
    cited = leading.group(1).split(",") if leading else []
    loose = []
    for locant, n in sorted(bonding.items(), key=lambda it: _locant_key(it[0])):
        if locant in cited:
            cited[cited.index(locant)] += f"\u03bb{n}"
        else:
            loose.append(f"{locant}\u03bb{n}")
    body = (",".join(cited) + "-" + leading.group(2)) if cited else stem
    return body, (",".join(loose) + "-") if loose else ""


def _chalcogen_heterones(mol, skeleton_atoms):
    """Ring sulfur, selenium or tellurium atoms bearing doubly bonded oxygen: the λ4/λ6 atoms of heterone names (P-64.4.2)."""
    ring_info = mol.GetRingInfo()
    return {
        a
        for a in skeleton_atoms
        if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52)
        and ring_info.NumAtomRings(a) == 1
        and not mol.GetAtomWithIdx(a).GetFormalCharge()
        and any(
            n.GetIdx() not in skeleton_atoms
            and n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(a, n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in mol.GetAtomWithIdx(a).GetNeighbors()
        )
    }


def _fused_mancude(mol, skeleton_atoms):
    oxo_all = _exocyclic_oxo(mol, skeleton_atoms)
    heterones = _chalcogen_heterones(mol, skeleton_atoms) & SUFFIX_ATOMS.get()
    oxo_all = oxo_all | heterones
    ring_info = mol.GetRingInfo()
    halogens = _ring_halogens(mol, skeleton_atoms)
    sp3 = {
        a
        for a in _sp3_ring_atoms(mol, skeleton_atoms) | oxo_all
        if mol.GetAtomWithIdx(a).GetSymbol() not in _NO_DOUBLE_BOND or a in halogens
    }
    sp3 |= {a for a in ANION_SUFFIX.get() if a in skeleton_atoms and mol.GetAtomWithIdx(a).GetAtomicNum() == 6}
    sp3 |= heterones
    fusion_hetero = {a for a in skeleton_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and ring_info.NumAtomRings(a) > 1}
    suffix_atoms = SUFFIX_ATOMS.get() & set(skeleton_atoms)
    oxo_suffix = (oxo_all & suffix_atoms) | {a for a in suffix_atoms & sp3 if ring_info.NumAtomRings(a) > 1}
    bare, old_of, parent = _named_mancude(mol, skeleton_atoms, sp3)
    match = re.match(r"^(\d+H(?:,\d+H)*)-(.*)$", parent)
    ih_count = len(match.group(1).split(",")) if match else 0
    stem = match.group(2) if match else parent
    if halogens:
        stem = re.sub(r"(?<=\d)\u03bb\d+", "", stem)
    numberings = fusion_system_numberings(bare, stem if stem in EXCEPTIONS else None)
    replacement = _replacement_numberings(bare)
    if replacement:
        numberings = replacement
    if not numberings:
        raise UnsupportedStructure("this fused skeleton has no supported peripheral numbering as a diyl yet")
    adj, can_hold = _ring_graph(mol, skeleton_atoms)
    can_hold -= fusion_hetero
    if heterones or halogens:
        can_hold |= heterones
        lambda_parent = nx.Graph()
        lambda_parent.add_nodes_from(can_hold)
        lambda_parent.add_edges_from((a, b) for a in can_hold for b in adj[a] if b in can_hold)
        ih_count = len(can_hold) - 2 * len(nx.max_weight_matching(lambda_parent, maxcardinality=True))
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
                    and not mol.GetAtomWithIdx(a).HasProp("_ring_cation_centre")
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
        lam = {a: n for a in position_of if (n := nonstandard_bonding(mol.GetAtomWithIdx(a)))}
        lam.update({a: int(mol.GetAtomWithIdx(a).GetTotalValence()) for a in heterones if a in position_of})
        marked_stem = stem
        if lam:
            marked_stem, loose = _cite_lambda(stem, {str(position_of[a]): n for a, n in lam.items()})
            ih_text += loose
        delta = _delta_citation(mol, skeleton_atoms, position_of)
        delta_stem = ("Δ" + ",".join(t for _, t in delta) + "-" + marked_stem) if delta else marked_stem

        fully_saturated = bool(hydro) and can_hold <= saturated | set(accommodated)

        def text(locants, valence, substituted=frozenset(), suffix="yl", hydro=hydro, ih_text=ih_text, added=added, full=fully_saturated, stem=delta_stem):
            rest = ih_text + _tail_added(stem, locants, valence, suffix, added)
            hydro_text = multiplied_word(len(hydro), "hydro") if full else _hydro_text(hydro)
            return hydro_text + ("-" if hydro and rest[0].isdigit() else "") + rest

        lam_key = tuple(sorted((-n, position_of[a]) for a, n in lam.items()))
        numbering = Numbering(position_of, text, pre_key=(ih, lam_key), unsat_key=(tuple(p for p, _ in delta), added, hydro), ih=ih)
        numbering.hydro, numbering.added, numbering.fully_hydro, numbering.stem = tuple(hydro), tuple(added), fully_saturated, stem
        numbering.parent_stem, numbering.hydro_positions, numbering.added_positions = stem, tuple(hydro), tuple(added)
        out.append(numbering)
    if not out:
        raise UnsupportedStructure("no numbering of this partly hydrogenated fused system fits its hydro/indicated hydrogen")
    return out


def _bridged_numberings(mol, skeleton_atoms):
    """Numberings of a fused ring system with bridges (P-25.4): hydro prefixes and indicated hydrogen follow the
    maximum number of noncumulative double bonds left after the bridge bonds are allowed for (P-25.7.1.1)."""
    from ._fusion_bridged import bridged_parents

    parents = bridged_parents(mol, set(skeleton_atoms))
    oxo_all = _exocyclic_oxo(mol, skeleton_atoms)
    suffix_atoms = SUFFIX_ATOMS.get() & set(skeleton_atoms)
    out = []
    for parent in parents:
        capable = {a for a in parent.capable if mol.GetAtomWithIdx(a).GetSymbol() not in _NO_DOUBLE_BOND} - set(parent.consumed)
        capable = {a for a in capable if mol.GetAtomWithIdx(a).GetDegree() <= 3 or mol.GetAtomWithIdx(a).GetSymbol() == "C"}
        capable = {a for a in capable if not (mol.GetAtomWithIdx(a).GetAtomicNum() in _TRIVALENT_RING_HETERO and mol.GetAtomWithIdx(a).GetDegree() == 3)}
        adj = {a: {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in capable} for a in capable}
        sp3 = {
            a
            for a in capable
            if not any(
                mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble() != 1.0 and n in capable for n in adj[a]
            )
        }
        position_of = {a: _locant(str(p)) if not isinstance(p, int) else p for a, p in parent.position_of.items()}
        if any(position_of.get(a) is None for a in capable):
            continue
        graph_nx = nx.Graph()
        graph_nx.add_nodes_from(capable)
        graph_nx.add_edges_from((a, b) for a in adj for b in adj[a])
        ih_count = len(capable) - 2 * len(nx.max_weight_matching(graph_nx, maxcardinality=True)) - len(oxo_all & capable & suffix_atoms) * 0
        accommodated = oxo_all & suffix_atoms
        split = _split_hydrogen(position_of, adj, capable, sp3 | (oxo_all & capable), oxo_all & capable, accommodated & capable, ih_count)
        if split is None:
            continue
        ih, added, hydro = split
        ih_text = (",".join(f"{p}H" for p in ih) + "-") if ih else ""
        stem = parent.stem

        def text(locants, valence, substituted=frozenset(), suffix="yl", hydro=hydro, ih_text=ih_text, added=added, stem=stem):
            rest = ih_text + _tail_added(stem, locants, valence, suffix, added)
            hydro_text = _hydro_text(hydro)
            return hydro_text + ("-" if hydro and rest[0].isdigit() else "") + rest

        numbering = Numbering(position_of, text, pre_key=(ih,), unsat_key=(added, hydro), ih=ih)
        numbering.hydro, numbering.added, numbering.stem = tuple(hydro), tuple(added), stem
        out.append(numbering)
    if not out:
        raise UnsupportedStructure("no numbering of this bridged fused system fits its hydro and indicated hydrogen")
    return out


def bridged_ring_system_name(mol):
    """Name of a bare ring system made of a fused ring system and bridges (P-25.4), or None."""
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetRingInfo().NumRings() < 3:
        return None
    info = mol.GetRingInfo()
    if any(info.NumAtomRings(i) == 0 for i in range(mol.GetNumAtoms())):
        return None
    atoms = list(range(mol.GetNumAtoms()))
    if any(not bond.IsInRing() for bond in mol.GetBonds()):
        return None
    if not is_bridged_fusion_system(mol, atoms):
        return None
    try:
        numberings = _bridged_numberings(mol, atoms)
        best = min(numberings, key=lambda n: (n.pre_key, n.unsat_key))
        return best.text((), 0, frozenset(), "")
    except UnsupportedStructure:
        return None


def is_bridged_fusion_system(mol, skeleton_atoms):
    """P-52.2.4.4, P-52.2.5.2.1: a ring system of at least three rings is named as a bridged fused system when a fused ring
    system with two rings of five or more members remains after removing the bridges; `_fusion_bridged` finds it."""
    return mol.GetRingInfo().NumRings() >= 3


def is_hydro_fusion_system(mol, skeleton_atoms):
    """A non-aromatic ortho- or peri-fused ring system with at least two rings of five or more members
    is named as a hydro derivative of its mancude fusion parent (P-31.1.4.2.4), not by von Baeyer."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if set(r) <= set(skeleton_atoms)]
    if len(rings) < 2 or sum(len(r) >= 5 for r in rings) < 2:
        return False
    if any(a.GetIsAromatic() for a in map(mol.GetAtomWithIdx, skeleton_atoms)):
        return False
    pairs = [(i, j) for i in range(len(rings)) for j in range(i + 1, len(rings)) if rings[i] & rings[j]]
    if any(len(rings[i] & rings[j]) != 2 for i, j in pairs):
        return False
    reached, frontier = {0}, [0]
    while frontier:
        k = frontier.pop()
        for i, j in pairs:
            other = j if i == k else i if j == k else None
            if other is not None and other not in reached:
                reached.add(other)
                frontier.append(other)
    # ortho-fused (each atom in at most two rings) or peri-fused (an interior atom shared by three rings)
    return len(reached) == len(rings) and not any(sum(a in r for r in rings) > 3 for a in skeleton_atoms)


def _plain_fused(mol, skeleton_atoms):
    from ._fused_numbering import FusedSystem
    from ._fusion_components import skeleton

    order = sorted(skeleton_atoms)
    index = {a: i for i, a in enumerate(order)}
    sk = skeleton(
        [mol.GetAtomWithIdx(a).GetSymbol() for a in order],
        [(index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()]) for b in mol.GetBonds() if b.GetBeginAtomIdx() in index and b.GetEndAtomIdx() in index],
    )
    try:
        FusedSystem(sk)
    except UnsupportedStructure:
        return False
    return True


_APPENDIX3_MAPPED = {}


def _is_appendix3_system(mol, atoms):
    """A ring system that is, or lies in, an Appendix 3 retained parent (P-101): the parent is numbered as a whole, so
    the ring system gets no fusion numbering of its own."""
    from ._appendix3_skeletons import _MIN_SIZE, _best_skeleton

    if len(atoms) < _MIN_SIZE and mol.GetNumAtoms() < _MIN_SIZE:
        return False
    rings = sum(1 for r in mol.GetRingInfo().AtomRings() if set(r) <= set(atoms))
    key = Chem.MolToSmiles(mol)
    if key not in _APPENDIX3_MAPPED:
        try:
            best = _best_skeleton(mol, False)
        except Exception:
            best = None
        _APPENDIX3_MAPPED[key] = set(best[1].mapping.values()) if best else set()
    mapped = _APPENDIX3_MAPPED[key]
    if mapped and set(atoms) <= mapped and (mapped - set(atoms)) and rings >= 2:
        return True
    if len(atoms) < _MIN_SIZE or rings < 3:
        return False
    editable = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    fragment = editable.GetMol()
    try:
        sanitize_probe(fragment)
        return _best_skeleton(fragment, False) is not None
    except Exception:
        return False
def system_numberings(mol, graph, rings, skeleton_atoms):
    from ._polyspiro_union import has_spiro_union_shape, spiro_union_numberings

    if _is_appendix3_system(mol, skeleton_atoms):
        raise UnsupportedStructure("an Appendix 3 retained parent has its own numbering, not a fusion numbering (P-101)")
    if has_spiro_union_shape(mol, set(skeleton_atoms)):
        return spiro_union_numberings(mol, graph, skeleton_atoms)
    if is_bridged_fusion_system(mol, skeleton_atoms) and not _plain_fused(mol, skeleton_atoms):
        try:
            return _bridged_numberings(mol, skeleton_atoms)
        except UnsupportedStructure:
            pass
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in skeleton_atoms):
        return _fused_mancude(mol, skeleton_atoms)
    if is_hydro_fusion_system(mol, skeleton_atoms):
        try:
            return _fused_mancude(mol, skeleton_atoms)
        except UnsupportedStructure:
            pass
    if is_bridged_fusion_system(mol, skeleton_atoms):
        try:
            return _bridged_numberings(mol, skeleton_atoms)
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
