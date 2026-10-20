"""Glycosides, glycosyl halides, glycosylamines, disaccharides and oligosaccharides (P-102.4, P-102.5.6.2, P-102.6.1, P-102.7).

Each pyranose or furanose ring with its carbon tail is a unit. Units are cut at their bridging oxygens, named as plain
monosaccharides, and joined by glycosyl names with linkage locants; an anomeric aglycone, or an anomeric-anomeric bond,
turns the parent into a glycoside (ending 'ide'). A unit with several donors or groups on its oxygens is named as a
substituted monosaccharide with the donors as O-glycosyl prefixes.
"""

from rdkit import Chem

from ._carbohydrate import ring_sugar_name_and_order
from ._cited_group import cited_group, subtree
from ._common import InternalInvariantError, UnsupportedStructure, adjacency

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
    """({bridging oxygen: (unit, carbon, other unit, other carbon)}, (unit, oxygen, aglycone root, atoms) or None, [(unit,
    carbon, oxygen, root)] of the groups on the other oxygens of the units)."""
    owner = {c: u for u in units for c in u.carbons}
    bridges, aglycone, held = {}, None, []
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
                    held.append((u, c, o, other))
    return bridges, aglycone, held


def _unit_molecule(mol, graph, unit, groups, children):
    """(molecule, forced names) of `unit` with the groups on its oxygens (`groups`: (oxygen, root)) and a stand-in carbon
    for each glycosyl donor (`children`: (oxygen, text, compound)), the anomeric group turned into hydroxy."""
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    keep = set(unit.atoms)
    for oxygen, root in groups:
        keep |= subtree(graph, root, oxygen)
    forced = {}
    for k, (oxygen, text, compound) in enumerate(children):
        stand_in = editable.AddAtom(Chem.Atom(6))
        editable.GetAtomWithIdx(stand_in).SetIntProp("_orig", -1 - k)
        editable.AddBond(oxygen, stand_in, Chem.BondType.SINGLE)
        forced[-1 - k] = (text, compound)
    for index in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(index)
    for atom in editable.GetAtoms():
        if atom.GetIntProp("_orig") == unit.exo and unit.exo_kind != "oxygen":
            atom.SetAtomicNum(8)
        if atom.GetIntProp("_orig") not in unit.carbons and atom.GetIntProp("_orig") in keep:
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
    sub = editable.GetMol()
    Chem.SanitizeMol(sub)
    return sub, {a.GetIntProp("_orig"): a.GetIdx() for a in sub.GetAtoms()}, forced


def _substituted_name(mol, graph, unit, groups, children, ester_prefixes):
    from ._substituents import FORCED_BRANCH_NAMES
    from ._sugar_substituted import _substituted_sugar_name

    sub, index, forced = _unit_molecule(mol, graph, unit, groups, children)
    token = FORCED_BRANCH_NAMES.set((sub.GetNumAtoms(), {index[k]: v for k, v in forced.items()}))
    try:
        name = _substituted_sugar_name(sub, ester_prefixes=ester_prefixes)
    except InternalInvariantError:
        raise
    except (ValueError, RuntimeError, Chem.rdchem.MolSanitizeException):
        name = None
    finally:
        FORCED_BRANCH_NAMES.reset(token)
    if name is None:
        raise UnsupportedStructure("the groups on the oxygens of a sugar unit are not cited by any supported name")
    return name


def _reachable(unit, children):
    found, stack = {id(unit)}, [unit]
    while stack:
        for donor, _, _, _ in children.get(id(stack.pop()), []):
            found.add(id(donor))
            stack.append(donor)
    return found


def _name_tree(mol, graph, units):
    """(name, atoms of the groups the name cites besides the units, whether the name is a glycoside) of sugar units joined
    by glycosyl links (P-102.7); a link between two nonadjacent positions makes the unit beyond it a group of the unit
    on this side (P-102.6.2)."""
    bridges, aglycone, held = _bridges(mol, graph, units)
    if len(bridges) != len(units) - 1:
        raise UnsupportedStructure("the sugar units are not joined as a tree")
    parent_of, pair, ethers = {}, None, []
    for o, (u, c, v, n) in bridges.items():
        donor_u, donor_v = c == u.anomeric, n == v.anomeric
        if donor_u and donor_v:
            if pair is not None:
                raise UnsupportedStructure("more than one anomeric-anomeric bond")
            pair = (u, v)
        elif donor_u:
            parent_of[id(u)] = (v, n, u, o)
        elif donor_v:
            parent_of[id(v)] = (u, c, v, o)
        else:
            ethers.append((o, u, c, v, n))
    children = {}
    for parent, carbon, donor, oxygen in parent_of.values():
        children.setdefault(id(parent), []).append((donor, parent, carbon, oxygen))
    extra = set(aglycone[3]) if aglycone is not None else set()
    partner = word = None
    tops = [u for u in units if id(u) not in parent_of and u.exo_kind != "substituted amine" and mol.GetAtomWithIdx(u.exo).GetDegree() == 1]
    if aglycone is not None:
        root = aglycone[0]
        word = cited_group(mol, graph, aglycone[2], aglycone[1])[0]
    elif pair is not None:
        root, partner = sorted(pair, key=_parent_key)
    elif tops:
        functional = [u for u in tops if u.exo_kind != "oxygen"]
        if len(functional) > 1:
            raise UnsupportedStructure("more than one glycosyl halide or amine")
        root = functional[0] if functional else min(tops, key=_parent_key)
    else:
        raise UnsupportedStructure("no unique reducing unit")
    if id(root) in parent_of:
        raise UnsupportedStructure("the parent unit is also a glycosyl donor")
    main = _reachable(root, children) | (_reachable(partner, children) if partner is not None else set())
    groups, stood_in = {}, {}
    for u, c, o, other in held:
        groups.setdefault(id(u), []).append((o, other))
        extra |= subtree(graph, other, o)
    for o, u, c, v, n in ethers:
        host, link = (u, n) if id(u) in main else (v, c)
        stood_in.setdefault(id(host), []).append((o, *cited_group(mol, graph, link, o)))
        extra |= subtree(graph, link, o)

    def ending(unit, role, plain=(), links=()):
        if not plain and not links:
            return {"glycoside": _glycoside_ending, "glycose": _glycose, "glycosyl": _glycosyl}[role](unit)
        if unit.uronic:
            raise UnsupportedStructure("a uronic acid unit with groups on its oxygens")
        name = _substituted_name(mol, graph, unit, plain, links, role != "glycose")
        if role == "glycosyl":
            return name[: -len("e")] + "yl"
        return name[: -len("ose")] + "oside" if role == "glycoside" else name

    def text_of(unit, role):
        kids = children.get(id(unit), [])
        own = groups.get(id(unit), [])
        if len(kids) == 1 and not own and id(unit) not in stood_in:
            donor, parent, carbon, _ = kids[0]
            return f"{text_of(donor, 'glycosyl')}-({donor.locant(donor.anomeric)}→{parent.locant(carbon)})-{ending(unit, role)}"
        links = [(o, t, "-(" in t or " " in t) for donor, _, _, o in kids for t in [text_of(donor, "glycosyl")]]
        return ending(unit, role, own, links + stood_in.get(id(unit), []))

    glycosidic = aglycone is not None or partner is not None
    if root.exo_kind == "amine":
        return text_of(root, "glycosyl") + "amine", extra, False, bool(groups or stood_in)
    if root.exo_kind == "halide":
        return f"{text_of(root, 'glycosyl')} {_HALIDES[mol.GetAtomWithIdx(root.exo).GetAtomicNum()]}", extra, False, bool(groups or stood_in)
    text = text_of(root, "glycoside" if glycosidic else "glycose")
    if partner is not None:
        text = f"{text_of(partner, 'glycosyl')} {text}"
    if word is not None:
        text = f"{word} {text}"
    return text, extra, glycosidic, bool(groups or stood_in)


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
            if (
                not any(a.GetFormalCharge() or not a.IsInRing() or a.GetAtomicNum() != 6 for a in ring.GetAtoms())
                and Chem.GetMolFrags(ring) == (tuple(range(ring.GetNumAtoms())),)
                and len(set(Chem.CanonicalRankAtoms(ring, breakTies=False))) == 1
            ):
                from .core import smiles_to_iupac

                return f"({name}){smiles_to_iupac(Chem.MolToSmiles(ring))}" if rest else None
            return _glycosyl_on_parent(mol, link, sugar, name)
    return None


def _glycosyl_on_parent(mol, link, sugar, name):
    """The parent that carries a glycosyl group on the carbon `link` (P-102.6.1.4): the sugar is swapped for a
    placeholder substituent whose name is supplied while the rest is named."""
    from .core import _name_mol
    from ._substituents import FORCED_BRANCH_NAMES

    editable = Chem.RWMol(mol)
    placeholder = editable.AddAtom(Chem.Atom(6))
    editable.AddBond(link.GetIdx(), placeholder, Chem.BondType.SINGLE)
    for index in sorted(sugar, reverse=True):
        editable.RemoveAtom(index)
    reduced = editable.GetMol()
    try:
        Chem.SanitizeMol(reduced)
    except Exception:
        return None
    token = FORCED_BRANCH_NAMES.set((reduced.GetNumAtoms(), {placeholder - len(sugar): (name, True)}))
    try:
        found = _name_mol(reduced)
    except UnsupportedStructure:
        return None
    finally:
        FORCED_BRANCH_NAMES.reset(token)
    return found if name in found and "methyl" not in found.replace(name, "") else None


def glycoside_name(mol):
    from ._substituents import ISOTOPE_LABELS

    if any(a.GetIsotope() for a in mol.GetAtoms()) or ISOTOPE_LABELS.get() is not None:
        return None
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
        text, outside, glycosidic, grouped = _name_tree(mol, graph, units)
        if len(units) == 1 and (grouped or (not glycosidic and not units[0].uronic)):
            return None
        return text if covered | outside == set(range(mol.GetNumAtoms())) else None
    except UnsupportedStructure:
        return None


def has_glycoside_shape(mol) -> bool:
    return glycoside_name(mol) is not None


def name_glycoside(mol) -> str:
    return glycoside_name(mol)
