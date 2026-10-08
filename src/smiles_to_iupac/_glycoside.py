"""Glycosides, glycosyl halides, glycosylamines, disaccharides and oligosaccharides (P-102.4, P-102.5.6.2, P-102.6.1, P-102.7).

Each pyranose or furanose ring with its carbon tail is a unit. Units are cut at their bridging oxygens, named as plain
monosaccharides, and joined by glycosyl names with linkage locants; an anomeric aglycone, or an anomeric-anomeric bond,
turns the parent into a glycoside (ending 'ide').
"""

from rdkit import Chem

from ._carbohydrate import ring_sugar_name_and_order
from ._cited_group import cited_group, subtree
from ._common import UnsupportedStructure, adjacency

_HALIDES = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_AGLYCONE_ELEMENTS = {6, 8, 9, 17, 35, 53}
_MAX_CHAIN = 2
_URONIC = {"acid": ("uronic acid",), "amide": ("uronamide",)}


class _Unit:
    def __init__(self, atoms, carbons, anomeric, exo, exo_kind, name, order, uronic=False):
        self.uronic = uronic
        self.atoms = atoms
        self.carbons = carbons
        self.anomeric = anomeric
        self.exo = exo
        self.exo_kind = exo_kind
        self.name = name
        self.order = order

    def locant(self, carbon):
        return self.order.index(carbon) + 1

    @property
    def is_ketose(self):
        return self.order[1] == self.anomeric


def _chain_carbons(mol, ring_set, ring_carbons):
    chain, frontier = set(), list(ring_carbons)
    while frontier:
        for n in mol.GetAtomWithIdx(frontier.pop()).GetNeighbors():
            idx = n.GetIdx()
            if idx in ring_set or idx in chain or n.GetAtomicNum() != 6 or n.IsInRing():
                continue
            chain.add(idx)
            frontier.append(idx)
    return chain


def _anomeric_carbon(mol, ring_set, ring_oxygen):
    candidates = [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(ring_oxygen).GetNeighbors()
        if any(m.GetIdx() not in ring_set and m.GetAtomicNum() != 6 for m in n.GetNeighbors())
    ]
    return candidates[0] if len(candidates) == 1 else None


def _exo_kind(atom):
    if atom.GetAtomicNum() == 8:
        return "oxygen"
    if atom.GetAtomicNum() == 7 and not atom.GetFormalCharge() and not atom.IsInRing():
        if atom.GetDegree() == 1 and atom.GetTotalNumHs() == 2:
            return "amine"
        if atom.GetDegree() == 2 and atom.GetTotalNumHs() == 1:
            return "substituted amine"
    if atom.GetAtomicNum() in _HALIDES and atom.GetDegree() == 1:
        return "halide"
    return None


def _capped_unit(mol, ring):
    ring_set = set(ring)
    ring_oxygens = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]
    if len(ring_oxygens) != 1:
        return None
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() not in (6, 8) or mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return None
    ring_carbons = ring_set - set(ring_oxygens)
    chain = _chain_carbons(mol, ring_set, ring_carbons)
    if len(chain) > _MAX_CHAIN:
        return None
    carbons = ring_carbons | chain
    atoms = ring_set | chain
    for c in carbons:
        for n in mol.GetAtomWithIdx(c).GetNeighbors():
            if n.GetIdx() in atoms:
                continue
            if n.GetAtomicNum() not in (7, 8, 9, 17, 35, 53):
                return None
            atoms.add(n.GetIdx())
    anomeric = _anomeric_carbon(mol, ring_set, ring_oxygens[0])
    if anomeric is None:
        return None
    exo = next(m for m in mol.GetAtomWithIdx(anomeric).GetNeighbors() if m.GetIdx() not in ring_set and m.GetAtomicNum() != 6)
    kind = _exo_kind(exo)
    if kind is None:
        return None
    carbonyls = [
        n.GetIdx()
        for c in chain
        for n in mol.GetAtomWithIdx(c).GetNeighbors()
        if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(c, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if len(carbonyls) > 1:
        return None
    amide_nitrogens = [
        n.GetIdx()
        for c in chain
        for n in mol.GetAtomWithIdx(c).GetNeighbors()
        if carbonyls and n.GetAtomicNum() == 7 and n.GetDegree() == 1 and n.GetTotalNumHs() == 2 and not n.GetFormalCharge()
    ]
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    for index in sorted(set(range(mol.GetNumAtoms())) - atoms | set(carbonyls), reverse=True):
        editable.RemoveAtom(index)
    for atom in editable.GetAtoms():
        origin = atom.GetIntProp("_orig")
        if (origin == exo.GetIdx() and kind != "oxygen") or origin in amide_nitrogens:
            atom.SetAtomicNum(8)
        if origin not in ring_set and origin not in carbons:
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
    sub = editable.GetMol()
    try:
        Chem.SanitizeMol(sub)
    except Exception:
        return None
    found = ring_sugar_name_and_order(sub)
    if found is None:
        return None
    name, order = found
    origin = {a.GetIdx(): a.GetIntProp("_orig") for a in sub.GetAtoms()}
    unit = _Unit(atoms, carbons, anomeric, exo.GetIdx(), kind, name, tuple(origin[i] for i in order), ("amide" if amide_nitrogens else "acid") if carbonyls else False)
    return None if unit.uronic and unit.is_ketose else unit


def _units(mol):
    units, taken = [], set()
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6) or any(a in taken for a in ring):
            continue
        unit = _capped_unit(mol, ring)
        if unit is not None:
            units.append(unit)
            taken |= set(ring)
    return units


def _aglycone_ok(mol, graph, root, oxygen):
    if mol.GetAtomWithIdx(root).GetAtomicNum() != 6:
        return False
    for a in subtree(graph, root, oxygen):
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() not in _AGLYCONE_ELEMENTS or atom.GetFormalCharge() or atom.GetIsotope():
            return False
        if atom.GetAtomicNum() == 8 and (atom.GetDegree() > 2 or any(b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds())):
            return False
    return True


def _parent_key(unit):
    """P-102.4 (a)-(c): aldose before ketose, then the alphabetically first stem, D before L, alpha before beta."""
    anomeric = unit.name[0] if unit.name[0] in "αβ" else ""
    stem = unit.name[len(anomeric) + 1 if anomeric else 0 :]
    return (unit.is_ketose, stem[2:], stem[0] != "D", anomeric != "α")


def _glycose(unit):
    return unit.name[: -len("ose")] + _URONIC[unit.uronic][0] if unit.uronic else unit.name


def _glycosyl(unit):
    return unit.name[: -len("ose")] + "osyl" + _URONIC[unit.uronic][0] if unit.uronic else unit.name[: -len("e")] + "yl"


def _glycoside_ending(unit):
    return unit.name[: -len("ose")] + "osid" + _URONIC[unit.uronic][0] if unit.uronic else unit.name[: -len("e")] + "ide"


def _bridges(mol, graph, units):
    """({bridging oxygen: (unit, carbon, other unit, other carbon)}, (unit, oxygen, aglycone root, atoms) or None)."""
    owner = {c: u for u in units for c in u.carbons}
    bridges, aglycone = {}, None
    for u in units:
        for c in u.carbons:
            for n in mol.GetAtomWithIdx(c).GetNeighbors():
                o = n.GetIdx()
                if o in u.carbons or n.IsInRing() or o not in u.atoms:
                    continue
                others = [m.GetIdx() for m in n.GetNeighbors() if m.GetIdx() != c]
                if not others:
                    continue
                if n.GetAtomicNum() != 8 or len(others) != 1:
                    raise UnsupportedStructure("a substituted exocyclic atom outside the supported links")
                other = others[0]
                if other in owner and owner[other] is not u:
                    bridges[o] = (u, c, owner[other], other)
                elif c == u.anomeric and aglycone is None and _aglycone_ok(mol, graph, other, o):
                    aglycone = (u, o, other, subtree(graph, other, o))
                else:
                    raise UnsupportedStructure("a substituent on an oxygen of a sugar unit")
    return bridges, aglycone


def _name_tree(mol, graph, units):
    bridges, aglycone = _bridges(mol, graph, units)
    if len(bridges) != len(units) - 1:
        raise UnsupportedStructure("the sugar units are not joined as a tree")
    if len(units) > 1 and any(u.exo_kind != "oxygen" for u in units):
        raise UnsupportedStructure("a glycosyl halide or amine of an oligosaccharide")
    parent_of, pair = {}, None
    for u, c, v, n in bridges.values():
        donor_u, donor_v = c == u.anomeric, n == v.anomeric
        if donor_u and donor_v:
            if pair is not None:
                raise UnsupportedStructure("more than one anomeric-anomeric bond")
            pair = (u, v)
        elif donor_u:
            parent_of[id(u)] = (v, n, u)
        elif donor_v:
            parent_of[id(v)] = (u, c, v)
        else:
            raise UnsupportedStructure("a link between two non-anomeric positions")
    children = {}
    for parent, carbon, donor in parent_of.values():
        children.setdefault(id(parent), []).append((donor, parent, carbon))
    if any(len(below) > 1 for below in children.values()):
        raise UnsupportedStructure("branched oligosaccharide")
    partner = word = None
    if aglycone is not None:
        root = aglycone[0]
        word = cited_group(mol, graph, aglycone[2], aglycone[1])[0]
    elif pair is not None:
        root, partner = sorted(pair, key=_parent_key)
    else:
        free = [u for u in units if id(u) not in parent_of and len(mol.GetAtomWithIdx(u.exo).GetNeighbors()) == 1]
        if len(free) != 1:
            raise UnsupportedStructure("no unique reducing unit")
        root = free[0]
    if id(root) in parent_of:
        raise UnsupportedStructure("the parent unit is also a glycosyl donor")

    def chain(unit, role):
        ending = {"glycoside": _glycoside_ending, "glycose": _glycose, "glycosyl": _glycosyl}[role](unit)
        below = children.get(id(unit))
        if not below:
            return ending
        donor, parent, carbon = below[0]
        return f"{chain(donor, 'glycosyl')}-({donor.locant(donor.anomeric)}→{parent.locant(carbon)})-{ending}"

    text = chain(root, "glycose" if aglycone is None and partner is None else "glycoside")
    if partner is not None:
        text = f"{chain(partner, 'glycosyl')} {text}"
    if word is not None:
        text = f"{word} {text}"
    return text, (aglycone[3] if aglycone is not None else set())


def _c_glycosyl_ring(mol, graph):
    """'(glycosyl)ring' for a plain sugar bonded through its anomeric carbon to a symmetric unsubstituted ring (P-102.6.1.4)."""
    from ._glycosyl import glycosyl_group

    for atom in mol.GetAtoms():
        if not atom.IsInRing() or atom.GetAtomicNum() != 6:
            continue
        for link in atom.GetNeighbors():
            if link.IsInRing() and link.GetIdx() in subtree(graph, atom.GetIdx(), link.GetIdx()):
                continue
            found = glycosyl_group(mol, graph, atom.GetIdx(), link.GetIdx())
            if found is None or link.GetAtomicNum() != 6:
                continue
            name, sugar = found
            rest = set(range(mol.GetNumAtoms())) - sugar
            editable = Chem.RWMol(mol)
            for index in sorted(sugar, reverse=True):
                editable.RemoveAtom(index)
            ring = editable.GetMol()
            Chem.SanitizeMol(ring)
            if any(a.GetFormalCharge() or not a.IsInRing() or a.GetAtomicNum() != 6 for a in ring.GetAtoms()):
                return None
            if Chem.GetMolFrags(ring) != (tuple(range(ring.GetNumAtoms())),) or len(set(Chem.CanonicalRankAtoms(ring, breakTies=False))) != 1:
                return None
            from .core import smiles_to_iupac

            return f"({name}){smiles_to_iupac(Chem.MolToSmiles(ring))}" if rest else None
    return None


def glycoside_name(mol):
    if not any(len(r) in (5, 6) and sum(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in r) == 1 for r in mol.GetRingInfo().AtomRings()):
        return None
    try:
        units = _units(mol)
        if not units:
            return _c_glycosyl_ring(mol, adjacency(mol))
        graph = adjacency(mol)
        covered = set().union(*(u.atoms for u in units))
        if len(units) == 1 and units[0].exo_kind != "oxygen":
            unit = units[0]
            if unit.uronic:
                return None
            if covered != set(range(mol.GetNumAtoms())) and unit.exo_kind != "substituted amine":
                return None
            if unit.exo_kind == "amine":
                return f"{_glycosyl(unit)}amine"
            if unit.exo_kind == "substituted amine":
                root = next(n.GetIdx() for n in mol.GetAtomWithIdx(unit.exo).GetNeighbors() if n.GetIdx() not in unit.carbons)
                if not _aglycone_ok(mol, graph, root, unit.exo):
                    return None
                covered |= subtree(graph, root, unit.exo)
                if covered != set(range(mol.GetNumAtoms())):
                    return None
                return f"N-{cited_group(mol, graph, root, unit.exo)[0]}-{_glycosyl(unit)}amine"
            return f"{_glycosyl(unit)} {_HALIDES[mol.GetAtomWithIdx(unit.exo).GetAtomicNum()]}"
        text, aglycone_atoms = _name_tree(mol, graph, units)
        if len(units) == 1 and not aglycone_atoms and not units[0].uronic:
            return None
        return text if covered | aglycone_atoms == set(range(mol.GetNumAtoms())) else None
    except UnsupportedStructure:
        return None


def has_glycoside_shape(mol) -> bool:
    return glycoside_name(mol) is not None


def name_glycoside(mol) -> str:
    return glycoside_name(mol)
