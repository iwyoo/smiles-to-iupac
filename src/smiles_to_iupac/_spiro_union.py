"""Components of spiro unions of rings (P-24.3, P-24.5, P-24.8): monocycles, fused systems named by fusion nomenclature,
and polycyclic von Baeyer systems, each with the numberings its own name allows. Heteroatoms of von Baeyer components, and
of monocycles without a Hantzsch-Widman name, are cited as 'a' prefixes (P-24.5.2, P-24.3.4); the spiro atoms, hydrogen
and locant citation are decided in `_polyspiro_union.py`.
"""

import re
from itertools import combinations

import networkx as nx
from rdkit import Chem

from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, ring_cycle
from ._fused_numbering import HETERO_RANK
from ._fusion_name import Context, fusion_name, system_numbering_options
from ._hetero_monocyclic import has_hetero_monocyclic_name, name_hetero_monocyclic
from ._numerals import alkane_name
from ._ring_diyl_numbering import _MANCUDE_RETAINED, _RANK, _hantzsch_widman, _saturated_name, _with_hetero_locants
from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates
from ._von_baeyer_heteroatom import _replacement_multiplied_word

_PRIME = "′"
_STANDARD = {5: 3, 7: 3, 8: 2, 15: 3, 16: 2, 33: 3, 34: 2, 51: 3, 52: 2, 83: 3}
_GROUP_14 = {6, 14, 32, 50, 82}
_HALOGENS = {9, 17, 35, 53}
_A_PREFIX = {
    8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza", 15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma",
    14: "sila", 32: "germa", 50: "stanna", 82: "plumba", 5: "bora",
}
_MAX_HANTZSCH_WIDMAN_SIZE = 10
_MIN_FUSION_RING = 5


class SpiroLocant:
    """A locant of a spiro system, ordered by number, then primes, then fusion letters: 4 < 4a < 4′ < 4′a < 5."""

    def __init__(self, text):
        match = re.fullmatch(r"(\d+)(′*)([a-z]*)", text)
        self.text = text
        self.key = (int(match.group(1)), len(match.group(2)), match.group(3))

    def __eq__(self, other):
        return isinstance(other, SpiroLocant) and self.key == other.key

    def __lt__(self, other):
        return self.key < other.key

    def __le__(self, other):
        return self.key <= other.key

    def __gt__(self, other):
        return self.key > other.key

    def __ge__(self, other):
        return self.key >= other.key

    def __hash__(self):
        return hash(self.key)

    def __str__(self):
        return self.text

    __repr__ = __str__

    def __format__(self, spec):
        return self.text


def _lk(text):
    match = re.match(r"(\d+)(′*)([a-z]*)", text)
    return int(match.group(1)), len(match.group(2)), match.group(3)


def _primed(locant, count):
    match = re.match(r"(\d+)(.*)", locant)
    return match.group(1) + _PRIME * count + match.group(2)


def _components(mol):
    info = mol.GetRingInfo()
    rings = [set(r) for r in info.AtomRings()]
    bond_rings = [set(r) for r in info.BondRings()]
    parent = list(range(len(rings)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in combinations(range(len(rings)), 2):
        if bond_rings[i] & bond_rings[j]:
            parent[find(i)] = find(j)
    groups = {}
    for i in range(len(rings)):
        groups.setdefault(find(i), []).append(i)
    result = []
    for g in groups.values():
        atoms = set().union(*(rings[i] for i in g))
        bonds = set().union(*(bond_rings[i] for i in g))
        result.append(
            {
                "atoms": atoms,
                "bonds": bonds,
                "rings": len(bonds) - len(atoms) + 1,
                "von_baeyer": len(g) > 1
                and (
                    any(len(bond_rings[i] & bond_rings[j]) > 1 for i, j in combinations(g, 2))
                    or sum(len(rings[i]) >= _MIN_FUSION_RING for i in g) < 2
                ),
            }
        )
    return result


def _skeleton(mol, comp):
    atoms = sorted(comp["atoms"])
    index = {a: i for i, a in enumerate(atoms)}
    sub = Chem.RWMol()
    for a in atoms:
        sub.AddAtom(Chem.Atom(mol.GetAtomWithIdx(a).GetAtomicNum()))
    for b in comp["bonds"]:
        bond = mol.GetBondWithIdx(b)
        sub.AddBond(index[bond.GetBeginAtomIdx()], index[bond.GetEndAtomIdx()], Chem.BondType.SINGLE)
    result = sub.GetMol()
    Chem.SanitizeMol(result)
    Chem.GetSymmSSSR(result)
    return result, atoms


def _bracket_locants(name):
    match = re.match(r"(\d+(?:,\d+)*)-(.+)", name)
    return f"[{match.group(1)}]{match.group(2)}" if match else name


def _monocycle(sub, atoms, mol, force_replacement):
    double = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in sub.GetBonds()
        if mol.GetBondBetweenAtoms(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()]).GetBondTypeAsDouble() == 2.0
    }
    if any(mol.GetBondBetweenAtoms(atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()]).GetBondTypeAsDouble() not in (1.0, 2.0) for b in sub.GetBonds()):
        raise UnsupportedStructure("this bond order is not supported in a spiro component")
    hetero = any(a.GetAtomicNum() != 6 for a in sub.GetAtoms())
    replacement = hetero and (
        force_replacement or (sub.GetNumAtoms() > _MAX_HANTZSCH_WIDMAN_SIZE and not has_hetero_monocyclic_name(sub))
    )
    if replacement and double:
        raise UnsupportedStructure("an unsaturated monocyclic spiro component named by replacement is not supported yet")
    mancude = hetero and not replacement and bool(double)
    if mancude:
        name = None
    elif hetero and not replacement:
        name = name_hetero_monocyclic(sub) if has_hetero_monocyclic_name(sub) else None
    else:
        name = "cyclo" + alkane_name(sub.GetNumAtoms())
    cycle = ring_cycle(adjacency(sub), list(range(sub.GetNumAtoms())))
    size = len(cycle)
    options = []
    for start in range(size):
        for step in (1, -1):
            order = [cycle[(start + step * k) % size] for k in range(size)]
            locants = {a: k + 1 for k, a in enumerate(order)}
            symbols = [
                (HETERO_RANK.get(sub.GetAtomWithIdx(a).GetSymbol(), 99), locants[a])
                for a in order
                if sub.GetAtomWithIdx(a).GetAtomicNum() != 6 and not replacement
            ]
            ene = () if mancude else tuple(sorted(min(locants[a] for a in pair) for pair in double))
            key = (sorted(l for _, l in symbols), [l for _, l in sorted(symbols)])
            options.append((key, {a: str(l) for a, l in locants.items()}, ene))
    best = min(k for k, _, _ in options)
    kept = [(n, e) for k, n, e in options if k == best]
    if mancude:
        name = _mancude_name(sub, kept[0][0])
    elif name is None:
        name = _saturated_hetero_name(sub, kept[0][0])
    return name, kept, replacement, mancude


def _saturated_hetero_name(sub, numbering):
    by_locant = sorted(numbering, key=lambda a: int(numbering[a]))
    elements = tuple(sub.GetAtomWithIdx(a).GetSymbol() for a in by_locant)
    hetero = [(i + 1, _RANK.get(e, 99)) for i, e in enumerate(elements) if e != "C"]
    return _saturated_name(elements, hetero)


def _mancude_name(sub, numbering):
    by_locant = sorted(numbering, key=lambda a: int(numbering[a]))
    elements = tuple(sub.GetAtomWithIdx(a).GetSymbol() for a in by_locant)
    hetero = [(i + 1, _RANK.get(e, 99)) for i, e in enumerate(elements) if e != "C"]
    stem = _MANCUDE_RETAINED.get(elements)
    if stem is None:
        stem = _hantzsch_widman(elements, saturated=False)
        if stem is None:
            raise UnsupportedStructure("this unsaturated heteromonocycle has no supported mancude name")
        stem = _with_hetero_locants(stem, elements, hetero)
    return _bracket_locants(stem)


def _von_baeyer(sub):
    rings = sub.GetNumBonds() - sub.GetNumAtoms() + 1
    carbon = Chem.RWMol()
    for _ in range(sub.GetNumAtoms()):
        carbon.AddAtom(Chem.Atom(6))
    for bond in sub.GetBonds():
        carbon.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE)
    carbon = carbon.GetMol()
    if rings == 2:
        core = find_bicyclic_core(carbon)
        if core is None:
            raise UnsupportedStructure("this bridged spiro component has no von Baeyer name here")
        return bicyclic_parent_name(core), [{a: str(i + 1) for i, a in enumerate(order)} for order in iter_bicyclic_numberings(core)]
    core = find_polycyclic_core(carbon, rings) if 3 <= rings <= 6 else None
    if core is None:
        raise UnsupportedStructure("this bridged spiro component has no von Baeyer name here")
    candidates = list(iter_polycyclic_candidates(core, rings))
    best = min(outer for _, _, outer in candidates)
    names = {parent for _, parent, outer in candidates if outer == best}
    if len(names) != 1:
        raise UnsupportedStructure("this bridged spiro component has no von Baeyer name here")
    (name,) = names
    if name == "tricyclo[3.3.1.1^3,7]decane":
        raise UnsupportedStructure("adamantane is a retained name and is not supported in a spiro union yet")
    return name, [{a: str(i + 1) for i, a in enumerate(order)} for order, _, outer in candidates if outer == best]


def _component(mol, comp, force_replacement):
    sub, atoms = _skeleton(mol, comp)
    replacement = mancude = False
    if comp["rings"] == 1:
        name, numberings, replacement, mancude = _monocycle(sub, atoms, mol, force_replacement)
        if not replacement:
            name = _bracket_locants(name)
    elif comp["von_baeyer"]:
        replacement = True
        name, orders = _von_baeyer(sub)
        numberings = [(n, ()) for n in orders]
    elif force_replacement:
        raise UnsupportedStructure("a fused component beside a skeletal replacement spiro heteroatom is not supported yet")
    else:
        name, root = fusion_name(sub)
        numberings = [(n, ()) for n in system_numbering_options(Context(sub), name, root)]
        name = _bracket_locants(name).replace("'", _PRIME)
    return {
        "name": name,
        "numberings": [({atoms[i]: loc for i, loc in n.items()}, ene) for n, ene in numberings],
        "replacement": replacement,
        "fused": mancude or (comp["rings"] > 1 and not replacement),
    }


def _bonding_number(atom):
    z = atom.GetAtomicNum()
    standard = 4 if z in _GROUP_14 else _STANDARD.get(z)
    valence = atom.GetTotalValence()
    if standard is None or z not in _A_PREFIX or atom.GetFormalCharge() or valence < standard or (valence - standard) % 2:
        raise UnsupportedStructure("this skeletal replacement heteroatom is not supported")
    return valence if valence > standard else None


def _replacement_prefixes(mol, atoms, locant_of, with_lambda=True):
    by_element = {}
    for a in atoms:
        by_element.setdefault(mol.GetAtomWithIdx(a).GetAtomicNum(), []).append(a)
    groups = []
    for z in sorted(by_element, key=lambda z: HETERO_RANK[Chem.GetPeriodicTable().GetElementSymbol(z)]):
        cited = []
        for a in sorted(by_element[z], key=lambda a: _lk(locant_of[a])):
            lam = _bonding_number(mol.GetAtomWithIdx(a)) if with_lambda else None
            cited.append(f"{locant_of[a]}λ{lam}" if lam else locant_of[a])
        groups.append(f"{','.join(cited)}-{_replacement_multiplied_word(len(cited), _A_PREFIX[z])}")
    return "-".join(groups)


def _capable(mol, comp, spiro_atoms):
    result = set()
    for a in comp["atoms"]:
        if a in spiro_atoms:
            continue
        atom = mol.GetAtomWithIdx(a)
        ring_degree = sum(1 for n in atom.GetNeighbors() if mol.GetBondBetweenAtoms(a, n.GetIdx()).IsInRing())
        z = atom.GetAtomicNum()
        if z in _GROUP_14:
            capacity = 1 if ring_degree <= 3 else 0
        elif z in _STANDARD:
            if atom.GetTotalValence() - atom.GetFormalCharge() > _STANDARD[z]:
                raise UnsupportedStructure("a nonstandard bonding number away from the spiro atom is not supported yet")
            capacity = 1 if _STANDARD[z] - ring_degree == 1 else 0
        else:
            raise UnsupportedStructure("this ring atom is not supported in a spiro component yet")
        if capacity:
            result.add(a)
    return result


def _hydrogen_choices(mol, capable, polycyclic_bonds, monocyclic_bonds):
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    doubled = set()
    for bond in kekule.GetBonds():
        if bond.GetBondTypeAsDouble() == 2.0:
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bond.GetIdx() in monocyclic_bonds:
                continue
            if a not in capable or b not in capable or bond.GetIdx() not in polycyclic_bonds:
                raise UnsupportedStructure("a double bond outside the polycyclic spiro components is not supported yet")
            doubled |= {a, b}
        elif bond.GetBondTypeAsDouble() not in (1.0,):
            raise UnsupportedStructure("this bond order is not supported in a spiro system")
    saturated = capable - doubled
    graph = nx.Graph()
    graph.add_nodes_from(capable)
    for bond_idx in polycyclic_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in capable and b in capable:
            graph.add_edge(a, b)

    def unmatched(removed):
        sub = graph.subgraph(capable - set(removed))
        return len(sub) - 2 * len(nx.max_weight_matching(sub, maxcardinality=True))

    required = unmatched(())
    choices = [set(i) for i in combinations(sorted(saturated), required) if unmatched(i) == 0]
    if not choices:
        raise UnsupportedStructure("the indicated hydrogen of this spiro system cannot be placed")
    return saturated, choices


def _name_key(name):
    stripped = re.sub(r"-[\d,]+-", "", re.sub(r"^\[?\d+(?:,\d+)*\]?-?", "", name))
    lead = re.match(r"\[?(\d+(?:,\d+)*)\]?", name)
    locants = tuple(sorted(int(x) for x in lead.group(1).split(","))) if lead and re.match(r"\[?\d", name) else ()
    descriptor = re.search(r"cyclo\[([\d.^,]+)\]", name)
    numbers = tuple(int(x) for x in re.findall(r"\d+", descriptor.group(1))) if descriptor else ()
    while True:
        reduced = re.sub(r"\[[^\[\]]*\]", "", stripped)
        if reduced == stripped:
            break
        stripped = reduced
    return alpha_sort_key(stripped), numbers, locants, name


def _spiro_atom(mol, spiro):
    atom = mol.GetAtomWithIdx(spiro)
    z = atom.GetAtomicNum()
    charge = atom.GetFormalCharge()
    if z == 6:
        if charge or atom.GetTotalNumHs() or atom.GetDegree() != 4:
            raise UnsupportedStructure("this spiro carbon is not supported")
        return None, False
    if z in _GROUP_14 and not charge and atom.GetTotalValence() == 4:
        return None, False
    if z not in _STANDARD or charge not in (0, 1):
        raise UnsupportedStructure("this spiro atom is not supported")
    bonding = atom.GetTotalValence() + charge
    if bonding <= _STANDARD[z] or (bonding - _STANDARD[z]) % 2:
        raise UnsupportedStructure("this spiro heteroatom has no lambda bonding number")
    return bonding, bool(charge)
