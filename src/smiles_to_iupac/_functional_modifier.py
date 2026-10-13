"""Functional-modifier names of carbonyl derivatives (P-15.2.2, P-66.6.5, P-68.3.1): the aldehyde or ketone
followed by a class word, as in propanal oxime, acetone semicarbazone or butanal ethyl hemiacetal.
The derivative group is replaced by C=O, the parent carbonyl compound is named by the engine, and the
substituents of the leaving part are named as separate modifier words."""

import re
from collections import namedtuple

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._diester_ring_diyl import _system_of
from ._numerals import alkane_name, multiplying_prefix
from ._parent_locants import locant_key, parent_numberings
from ._ring_system_seniority import ring_seniority_key
from ._substituents import (
    format_mononuclear_prefixes,
    format_substituent_prefixes,
    multiplied_prefix,
    name_branch,
    wrap_marks,
)

_Group = namedtuple("_Group", "kind carbons leaving unit core stereo bond", defaults=(None, None))
_Acetal = namedtuple("_Acetal", "entries bridge elements term")

_CHALCOGEN_SYMBOLS = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_CHALCOGEN_WORDS = {8: "", 16: "thio", 34: "seleno", 52: "telluro"}
_RULES = {
    "oxime": "P-68.3.1.1.2",
    "hydrazone": "P-68.3.1.2.2",
    "azine": "P-68.3.1.2.3",
    "semicarbazone": "P-68.3.1.2.5",
    "acetal": "P-66.6.5.1",
    "chalcogen acetal": "P-66.6.5.3",
    "hemiacetal": "P-66.6.5.2",
    "chalcogen hemiacetal": "P-66.6.5.4",
}

_SENIOR_TO_CARBONYL = tuple(
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=[O,S,Se,Te])[O,S,Se,Te,N,F,Cl,Br,I]",
        "[CX3](=[NX2])[O,S,Se,Te,N]",
        "[CX2]#[NX1]",
        "[#6X2](=*)=*",
        "[#15,#16,#34,#52,#33,#51,#5,#14;$(*~[OX1])]-[#8,#7,#9,#17,#35,#53]",
    )
)
_NITRO = Chem.MolFromSmarts("[NX3+](=O)[O-]")
_ALDEHYDE = Chem.MolFromSmarts("[CX3;H1,H2]=[OX1]")
_KETONE = Chem.MolFromSmarts("[#6][CX3](=[OX1])[#6]")
# P-64.2.1.2: general nomenclature keeps these names for the unsubstituted structures
_GENERAL_NAMES = {
    Chem.CanonSmiles(smiles): name
    for smiles, name in (
        ("CC(C)=O", "acetone"),
        ("CC(=O)c1ccccc1", "acetophenone"),
        ("O=C(c1ccccc1)c1ccccc1", "benzophenone"),
        ("O=C1C=CC(=O)C=C1", "1,4-benzoquinone"),
    )
}


def functional_modifier_names(mol, name_parent):
    """[(name, rule)] for the molecule, or [] when it holds no derivative group; raises UnsupportedStructure
    when it holds one this module cannot name."""
    nitro = {i for match in mol.GetSubstructMatches(_NITRO) for i in match}
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        (a.GetFormalCharge() and a.GetIdx() not in nitro) or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return []
    context = _Context(mol)
    groups = _detect(context)
    if not groups:
        return []
    if len({group.kind for group in groups}) > 1:
        raise UnsupportedStructure("P-15.2.4: functional class names are not recommended for polyfunctional compounds")
    _require_plain_groups(mol, groups)
    if groups[0].kind == "azine":
        if len(groups) > 1:
            raise UnsupportedStructure("several azine groups")
        return [_azine_name(context, groups[0], name_parent)]
    parent, converted = _restore(mol, groups)
    principal = _principal_carbonyls(parent, converted)
    display, systematic = _parent_names(parent, name_parent)
    labels = {group.carbons[0]: group.stereo for group in groups if group.stereo}
    if labels:
        _require_citable_stereo(mol, groups)
    in_modifier = _locants_required(parent, converted, principal) or len({group.unit for group in groups}) > 1
    in_descriptor = bool(labels) and len(_descriptor_items(display)) + len(labels) > 1
    locants = _choose_locants(parent, groups, principal, systematic) if in_modifier or in_descriptor else None
    name = _with_descriptors(display, labels, locants)
    return [(f"{name} {_modifier(groups, locants if in_modifier else None)}", _rule(groups[0]))]


_DESCRIPTOR_ITEM = re.compile(r"(\d+[a-z]*'*)([RSEZ])")
_LEADING_DESCRIPTOR = re.compile(r"\(([^()]*)\)-")
_DESCRIPTOR_RANK = {"Z": 0, "E": 1}


def _descriptor_items(parent_name):
    match = _LEADING_DESCRIPTOR.match(parent_name)
    if match is None:
        return []
    items = match.group(1).split(",")
    if not all(_DESCRIPTOR_ITEM.fullmatch(item) for item in items):
        raise UnsupportedStructure("the stereodescriptors of the parent are not simple locant-descriptor pairs")
    return items


def _with_descriptors(parent_name, labels, locants):
    """P-93.4.2.1.3: the E/Z label of each C=N bond joins the descriptors of the parent, cited by the locant of its
    carbon; a lone descriptor of the whole name carries no locant."""
    if not labels:
        return parent_name
    items = _descriptor_items(parent_name)
    rest = parent_name[_LEADING_DESCRIPTOR.match(parent_name).end():] if items else parent_name
    if not items and len(labels) == 1:
        return f"({next(iter(labels.values()))})-{rest}"
    items += [f"{locants[carbon]}{label}" for carbon, label in labels.items()]
    items.sort(key=lambda item: locant_key(_DESCRIPTOR_ITEM.fullmatch(item).group(1)))
    return f"({','.join(items)})-{rest}"


def _require_citable_stereo(mol, groups):
    """Stereo elements of the parent are named with it; one inside a leaving part would be lost."""
    leaving = set().union(*(group.leaving for group in groups))
    own = {group.bond for group in groups}
    for element in Chem.FindPotentialStereo(mol):
        if element.specified != Chem.StereoSpecified.Specified:
            continue
        if element.type == Chem.StereoType.Bond_Double:
            if element.centeredOn in own:
                continue
            bond = mol.GetBondWithIdx(element.centeredOn)
            atoms = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
        else:
            atoms = {element.centeredOn}
        if atoms & leaving:
            raise UnsupportedStructure("a stereo element of the leaving part cannot be cited with the parent")


def _choose_locants(parent, groups, principal, systematic):
    """{derivative carbon: parent locant} of the numbering among the tied ones that cites the derivative groups
    lowest, then the substituent words first in alphanumerical order lowest (P-31.1.4.2.4, P-31.1.4.3.4)."""
    in_parent = {atom.GetIntProp("_origin"): atom.GetIdx() for atom in parent.GetAtoms() if atom.HasProp("_origin")}
    best = None
    for numbering in parent_numberings(parent, principal, systematic):
        if any(in_parent[c] not in numbering for group in groups for c in group.carbons):
            continue
        locants = {c: numbering[in_parent[c]] for group in groups for c in group.carbons}
        labelled = sorted((g for g in groups if g.stereo), key=lambda g: locant_key(locants[g.carbons[0]]))
        key = (
            sorted(locant_key(locants[c]) for c in locants),
            _citation_key(groups, locants),
            [_DESCRIPTOR_RANK[g.stereo] for g in labelled],
        )
        if best is None or key < best[0]:
            best = key, locants
    if best is None:
        raise UnsupportedStructure("the parent locants of the derivatized groups are not available")
    return best[1]


def _citation_key(groups, locants):
    parts, _ = _modifier_parts(groups)
    ordered = sorted(parts, key=lambda part: alpha_sort_key(part[2]))
    return [sorted(locant_key(locants[c]) for c in carbons) for _, carbons, _ in ordered]


class _Context:
    def __init__(self, mol):
        self.mol = mol
        self.graph = adjacency(mol)
        self.halogens = halogen_substituents(mol)
        self.aromatic = frozenset(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())

    def substituent(self, root, back):
        atom = self.mol.GetAtomWithIdx(root)
        if atom.GetAtomicNum() != 6 or any(
            bond.GetBondTypeAsDouble() == 2 and bond.GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16, 34, 52)
            for bond in atom.GetBonds()
        ):
            raise UnsupportedStructure("the substituent of the derivative group is not a hydrocarbyl group")
        return name_branch(self.graph, root, back, self.halogens, self.aromatic, mol=self.mol)


def _reach(graph, starts, blocked):
    seen = set(starts)
    stack = list(starts)
    while stack:
        for neighbour in graph[stack.pop()]:
            if neighbour not in seen and neighbour not in blocked:
                seen.add(neighbour)
                stack.append(neighbour)
    return seen


def _leaving(context, carbon, starts):
    """The part removed with the derivative group; it must be joined to the parent only through `carbon`."""
    fragment = _reach(context.graph, starts, {carbon})
    rest = [n for n in context.graph[carbon] if n not in starts]
    if fragment & _reach(context.graph, rest, {carbon}):
        raise UnsupportedStructure("a ring runs through the derivative group")
    return fragment


def _plain_carbon(atom, excluded):
    return (
        atom.GetAtomicNum() == 6
        and not atom.GetIsAromatic()
        and not atom.GetFormalCharge()
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors() if n.GetIdx() not in excluded)
    )


def _detect(context):
    mol = context.mol
    found = {}
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE or bond.GetIsAromatic():
            continue
        ends = (bond.GetBeginAtom(), bond.GetEndAtom())
        carbon = next((a for a in ends if a.GetAtomicNum() == 6), None)
        nitrogen = next((a for a in ends if a.GetAtomicNum() == 7), None)
        if carbon is None or nitrogen is None:
            continue
        group = _imine_derivative(context, carbon, nitrogen, bond)
        if group is not None:
            found[group.carbons] = group
    for atom in mol.GetAtoms():
        group = _acetal_group(context, atom)
        if group is not None:
            found[group.carbons] = group
    return list(found.values())


def _imine_derivative(context, carbon, nitrogen, bond):
    if not _plain_carbon(carbon, {nitrogen.GetIdx()}) or nitrogen.GetDegree() != 2 or nitrogen.GetFormalCharge():
        return None
    other = next(n for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon.GetIdx())
    if other.GetAtomicNum() not in (7, 8):
        return None
    c, n = carbon.GetIdx(), nitrogen.GetIdx()
    group = _oxime(context, c, n, other) if other.GetAtomicNum() == 8 else _nitrogen_derivative(context, c, n, other)
    label = _configuration(context.mol, bond)
    if label and group.kind == "azine":
        raise UnsupportedStructure("the configuration of an azine C=N bond is not cited")
    return group._replace(stereo=label, bond=bond.GetIdx())


def _configuration(mol, bond):
    """The CIP E or Z label of a stereodefined bond, with the lone pair of nitrogen ranking lowest (P-93.4.2.1.3)."""
    if bond.GetStereo() in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
        return None
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    return probe.GetBondWithIdx(bond.GetIdx()).GetProp("_CIPCode")


def _oxime(context, c, n, oxygen):
    if oxygen.GetFormalCharge() or oxygen.GetDegree() > 2:
        raise UnsupportedStructure("an oxime oxygen must be OH or O-substituted")
    substituent = None
    if oxygen.GetDegree() == 2:
        root = next(x for x in oxygen.GetNeighbors() if x.GetIdx() != n).GetIdx()
        substituent = context.substituent(root, oxygen.GetIdx())
    return _Group("oxime", (c,), _leaving(context, c, [n]), substituent, frozenset({c, n, oxygen.GetIdx()}))


def _nitrogen_derivative(context, c, n, other):
    if other.GetFormalCharge() or other.IsInRing() or context.mol.GetAtomWithIdx(n).IsInRing():
        raise UnsupportedStructure("a ring contains the hydrazone nitrogen")
    n2 = other.GetIdx()
    others = [x for x in other.GetNeighbors() if x.GetIdx() != n]
    double_to = [x for x in others if context.mol.GetBondBetweenAtoms(n2, x.GetIdx()).GetBondTypeAsDouble() == 2]
    if double_to:
        c2 = double_to[0]
        if len(others) != 1 or c2.GetAtomicNum() != 6 or not _plain_carbon(c2, {n2}):
            raise UnsupportedStructure("an azine nitrogen must be bonded to carbon only")
        pairs = ((c, n), (c2.GetIdx(), n2))
        return _Group("azine", tuple(sorted((c, c2.GetIdx()))), {n, n2}, pairs, frozenset({c, n, n2, c2.GetIdx()}))
    acyl = [x for x in others if _is_acyl(x)]
    if acyl:
        return _semicarbazone(context, c, n, other, acyl, others)
    entries = [context.substituent(x.GetIdx(), n2) for x in others]
    return _Group("hydrazone", (c,), _leaving(context, c, [n]), tuple(sorted(entries)), frozenset({c, n, n2}))


def _is_acyl(atom):
    return atom.GetAtomicNum() == 6 and any(
        b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(atom).GetAtomicNum() in _CHALCOGEN_SYMBOLS for b in atom.GetBonds()
    )


def _semicarbazone(context, c, n, nitrogen2, acyl, others):
    carbonyl = acyl[0]
    mol = context.mol
    chalcogen = [x for x in carbonyl.GetNeighbors() if mol.GetBondBetweenAtoms(carbonyl.GetIdx(), x.GetIdx()).GetBondTypeAsDouble() == 2]
    amide = [x for x in carbonyl.GetNeighbors() if x.GetIdx() != nitrogen2.GetIdx() and x.GetIdx() != chalcogen[0].GetIdx()]
    if len(acyl) != 1 or len(chalcogen) != 1 or len(amide) != 1 or amide[0].GetAtomicNum() != 7:
        raise UnsupportedStructure("an acylhydrazone other than a semicarbazone")
    n4 = amide[0]
    if n4.IsInRing() or n4.GetFormalCharge():
        raise UnsupportedStructure("a ring contains the semicarbazone nitrogen")
    entries = [(2,) + context.substituent(x.GetIdx(), nitrogen2.GetIdx()) for x in others if x.GetIdx() != carbonyl.GetIdx()]
    entries += [
        (4,) + context.substituent(x.GetIdx(), n4.GetIdx()) for x in n4.GetNeighbors() if x.GetIdx() != carbonyl.GetIdx()
    ]
    unit = (_CHALCOGEN_WORDS[chalcogen[0].GetAtomicNum()], tuple(sorted(entries)))
    core = frozenset({c, n, nitrogen2.GetIdx(), carbonyl.GetIdx(), chalcogen[0].GetIdx(), n4.GetIdx()})
    return _Group("semicarbazone", (c,), _leaving(context, c, [n]), unit, core)


def _acetal_group(context, atom):
    if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetFormalCharge():
        return None
    if any(bond.GetBondTypeAsDouble() != 1 for bond in atom.GetBonds()):
        return None
    ys = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in _CHALCOGEN_SYMBOLS]
    if len(ys) != 2 or any(n.GetTotalValence() != 2 or n.GetFormalCharge() for n in ys):
        return None
    if any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors() if n.GetIdx() not in {y.GetIdx() for y in ys}):
        return None
    c = atom.GetIdx()
    roots = [next((x.GetIdx() for x in y.GetNeighbors() if x.GetIdx() != c), None) for y in ys]
    if roots == [None, None]:
        return None
    y1, y2 = ys
    elements = tuple(sorted(y.GetAtomicNum() for y in ys))
    hemi = None in roots
    bridge = None
    entries = ()
    if not hemi and y2.GetIdx() in _reach(context.graph, [y1.GetIdx()], {c}):
        bridge = _bridge_name(context, c, y1.GetIdx(), y2.GetIdx())
        leaving = _leaving(context, c, [y1.GetIdx(), y2.GetIdx()])
    else:
        leaving = set()
        for y in ys:
            leaving |= _leaving(context, c, [y.GetIdx()])
        named = [None if r is None else context.substituent(r, y.GetIdx()) + (y.GetAtomicNum(),) for r, y in zip(roots, ys)]
        entries = tuple(sorted(e for e in named if e is not None))
    term = ("hemi" if hemi else "") + ("acetal" if atom.GetTotalNumHs() else "ketal")
    unit = _Acetal(entries, bridge, elements, term)
    kind = ("chalcogen " if elements != (8, 8) else "") + ("hemiacetal" if hemi else "acetal")
    return _Group(kind, (c,), leaving, unit, frozenset({c, y1.GetIdx(), y2.GetIdx()}))


def _chalcogen_class(elements):
    if elements == (8, 8):
        return ""
    if elements[0] == 8:
        return "mono" + _CHALCOGEN_WORDS[elements[1]]
    if elements[0] == elements[1]:
        return "di" + _CHALCOGEN_WORDS[elements[0]]
    return "".join(sorted(_CHALCOGEN_WORDS[e] for e in elements))


def _substituent_words(entries, elements):
    """Separate words for the O-, S-, Se- or Te-substituents, with the element symbol when the two chalcogens differ."""
    lettered = elements[0] != elements[1] and len({name for name, _, _ in entries}) == len(entries)
    if len(entries) == 2 and not lettered and entries[0][0] == entries[1][0]:
        return multiplied_prefix(2, entries[0][0], entries[0][1])
    ordered = sorted(entries, key=lambda e: alpha_sort_key(e[0]))
    words = []
    for name, compound, z in ordered:
        text = wrap_marks(name) if compound else name
        words.append(f"{_CHALCOGEN_SYMBOLS[z]}-{text}" if lettered else text)
    return " ".join(words)


def _bridge_name(context, c, y1, y2):
    """The divalent group joining the two chalcogens of a cyclic acetal, cited as a substituent group (P-66.6.5.1.2)."""
    mol = context.mol
    inner = _reach(context.graph, [y1], {c}) - {y1, y2}
    if not inner or any(
        mol.GetAtomWithIdx(a).GetAtomicNum() != 6 or mol.GetAtomWithIdx(a).GetDegree() != 2 or mol.GetAtomWithIdx(a).GetTotalNumHs() != 2
        for a in inner
    ):
        raise UnsupportedStructure("the group joining the acetal oxygens is not an unbranched alkylene chain")
    count = len(inner)
    if mol.GetAtomWithIdx(y1).GetDegree() != 2 or mol.GetAtomWithIdx(y2).GetDegree() != 2:
        raise UnsupportedStructure("the group joining the acetal oxygens is not an unbranched alkylene chain")
    return {1: "methylene", 2: "ethylene"}.get(count) or f"{alkane_name(count)}-1,{count}-diyl"


def _restore(mol, groups):
    """The molecule with every derivative group replaced by C=O, and the indices of the restored carbons."""
    rw = Chem.RWMol(mol)
    for c in (c for group in groups for c in group.carbons):
        rw.GetAtomWithIdx(c).SetIntProp("_origin", c)
        oxygen = rw.AddAtom(Chem.Atom(8))
        rw.AddBond(c, oxygen, Chem.BondType.DOUBLE)
    for index in sorted(set().union(*(group.leaving for group in groups)), reverse=True):
        rw.RemoveAtom(index)
    for atom in rw.GetAtoms():
        if atom.HasProp("_origin"):
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    parent = rw.GetMol()
    Chem.SanitizeMol(parent)
    return parent, _restored_atoms(parent)


def _locants_required(parent, converted, principal):
    """A single group on equivalent carbonyls, or every carbonyl, is cited without parent locants."""
    return not (len(converted) == len(principal) or (len(converted) == 1 and not _inequivalent(parent, principal)))


def _require_plain_groups(mol, groups):
    """Nothing outside the derivative groups may outrank the aldehyde or ketone they derive from (P-15.2.2: the
    modifier belongs to the principal characteristic group), and no isotope may be lost with the leaving part."""
    core = set().union(*(group.core for group in groups))
    for pattern in _SENIOR_TO_CARBONYL:
        if any(not set(match) <= core for match in mol.GetSubstructMatches(pattern)):
            raise UnsupportedStructure("a characteristic group senior to the aldehyde or ketone is present")
    if any(mol.GetAtomWithIdx(a).GetIsotope() for group in groups for a in group.leaving):
        raise UnsupportedStructure("an isotopically modified atom would be dropped with the derivative group")


def _principal_carbonyls(parent, converted):
    """The carbonyl carbons the parent name expresses by its suffix; the derivative groups must be among them."""
    principal = {m[0] for m in parent.GetSubstructMatches(_ALDEHYDE)} or {m[1] for m in parent.GetSubstructMatches(_KETONE)}
    if not converted <= principal:
        raise UnsupportedStructure("the derivative group is not the principal characteristic group of the parent")
    return principal


def _inequivalent(parent, principal):
    ranks = Chem.CanonicalRankAtoms(parent, breakTies=False)
    return len({ranks[a] for a in principal}) > 1


def _parent_names(parent, name_parent):
    """(name cited in the result, systematic name): general nomenclature keeps acetone and the like (P-64.2.1.2)."""
    smiles = Chem.MolToSmiles(parent)
    systematic = name_parent(smiles)
    return _GENERAL_NAMES.get(smiles) or systematic, systematic


def _rule(group):
    return _RULES[group.kind]


def _repeat(count, unit, simple):
    if count == 1:
        return unit
    if simple:
        return multiplying_prefix(count) + unit
    return multiplying_prefix(count, compound=True) + wrap_marks(unit)


_ACETAL_KINDS = ("acetal", "chalcogen acetal", "hemiacetal", "chalcogen hemiacetal")


def _modifier(groups, locants):
    """The modifier words after the parent; with `locants` ({carbon: parent locant}) each word is cited with the
    locants of the carbons it belongs to (P-66.6.5.1.1: '1,1-diethyl 4,4-dimethyl diketal')."""
    parts, tail = _modifier_parts(groups)
    parts.sort(key=lambda part: alpha_sort_key(part[2]))
    if locants is None:
        words = " ".join(text for text, _, _ in parts)
    else:
        if isinstance(groups[0].unit, _Acetal) and groups[0].unit.elements != (8, 8):
            raise UnsupportedStructure("the parent locants of chalcogen acetal carbons are not available")
        words = " ".join(_cited(text, carbons, locants) for text, carbons, _ in parts)
    return f"{words} {tail}" if tail else words


def _cited(text, carbons, locants):
    cited = ",".join(str(locants[c]) for c in sorted(carbons, key=lambda c: locant_key(locants[c])))
    return f"{cited}-{text if text.endswith(')') or not re.search(r'[- ]', text) else f'({text})'}"


def _modifier_parts(groups):
    """([(word, derivative carbons it applies to, word without multiplier)], class word shared by the words)."""
    kind = groups[0].kind
    carbons_of = {}
    for group in groups:
        carbons_of.setdefault(group.unit, []).extend(group.carbons)
    if kind not in _ACETAL_KINDS:
        return [_unit_part(*_unit_base(kind, unit), carbons) for unit, carbons in carbons_of.items()], ""
    if any(unit.bridge or unit.elements != (8, 8) for unit in carbons_of):
        return [_unit_part(_acetal_text(unit), False, carbons) for unit, carbons in carbons_of.items()], ""
    terms = {unit.term for unit in carbons_of}
    if len(terms) > 1:
        raise UnsupportedStructure("acetal and ketal carbons on one parent")
    by_name = {}
    for unit, carbons in carbons_of.items():
        for name, compound, _ in unit.entries:
            by_name.setdefault(name, (compound, []))[1].extend(carbons)
    parts = [(multiplied_prefix(len(c), name, compound), c, name) for name, (compound, c) in by_name.items()]
    term = terms.pop()
    return parts, (multiplying_prefix(len(groups)) + term if len(groups) > 1 else term)


def _acetal_text(unit):
    words = unit.bridge or _substituent_words(list(unit.entries), unit.elements)
    return f"{words} {_chalcogen_class(unit.elements)}{unit.term}"


def _unit_part(base, simple, carbons):
    return _repeat(len(carbons), base, simple), carbons, base


def _unit_base(kind, unit):
    """(text of one derivative group, whether a multiplier is added without enclosing marks)."""
    if kind == "oxime":
        if unit is None:
            return "oxime", True
        name, compound = unit
        return f"O-{wrap_marks(name) if compound else name}oxime", False
    if kind == "hydrazone":
        cluster = format_mononuclear_prefixes(list(unit)) if unit else ""
        return cluster + "hydrazone", not cluster
    chalcogen, entries = unit
    grouped = {}
    for locant, name, compound in entries:
        info = grouped.setdefault(name, {"locants": [], "compound": compound})
        info["locants"].append(locant)
    cluster = format_substituent_prefixes(grouped) if grouped else ""
    return cluster + chalcogen + "semicarbazone", not cluster


def _azine_name(context, group, name_parent):
    """P-68.3.1.2.3: a symmetrical azine takes the class word `azine`; an unsymmetrical one is the ylidenehydrazone
    of the senior aldehyde or ketone."""
    parent, _ = _restore(context.mol, [group])
    pieces = Chem.GetMolFrags(parent, asMols=True, sanitizeFrags=True)
    if len(pieces) != 2:
        raise UnsupportedStructure("a ring runs through the azine group")
    for piece in pieces:
        converted = _restored_atoms(piece)
        if _locants_required(piece, converted, _principal_carbonyls(piece, converted)):
            raise UnsupportedStructure("the parent locants of the azine carbons are not available")
    smiles = [Chem.MolToSmiles(piece) for piece in pieces]
    if smiles[0] == smiles[1]:
        return f"{_parent_names(pieces[0], name_parent)[0]} azine", _RULES["azine"]
    keys = [_component_key(piece) for piece in pieces]
    if keys[0] == keys[1]:
        raise UnsupportedStructure("the senior component of an unsymmetrical azine is not decided by class or ring")
    senior = keys.index(min(keys))
    carbon, nitrogen = next(
        pair for pair in group.unit if pair[0] == _origin(pieces[1 - senior])
    )
    ylidene, _ = name_branch(context.graph, carbon, nitrogen, context.halogens, context.aromatic, mol=context.mol)
    return f"{_parent_names(pieces[senior], name_parent)[0]} {ylidene}hydrazone", _RULES["azine"]


def _restored_atoms(parent):
    return {a.GetIdx() for a in parent.GetAtoms() if a.HasProp("_origin")}


def _origin(piece):
    return next(a.GetIntProp("_origin") for a in piece.GetAtoms() if a.HasProp("_origin"))


def _component_key(piece):
    """P-41 aldehyde before ketone, P-44.1.2.2 ring before chain, P-44.2 senior ring, P-44.3 longer chain."""
    carbon = next(a for a in piece.GetAtoms() if a.HasProp("_origin"))
    if not carbon.IsInRing() and piece.GetRingInfo().NumRings():
        raise UnsupportedStructure("a ring beside the carbonyl carbon: the parent hydride of the component is not decided here")
    aldehyde = 0 if carbon.GetTotalNumHs() else 1
    if carbon.IsInRing():
        return aldehyde, 0, ring_seniority_key(piece, _system_of(piece, carbon.GetIdx())[1])
    return aldehyde, 1, -_chain_length(piece, carbon.GetIdx())


def _chain_length(piece, carbon):
    def depth(index, back):
        arms = [depth(n.GetIdx(), index) for n in piece.GetAtomWithIdx(index).GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != back]
        return 1 + max(arms, default=0)

    arms = sorted((depth(n.GetIdx(), carbon) for n in piece.GetAtomWithIdx(carbon).GetNeighbors() if n.GetAtomicNum() == 6), reverse=True)
    return 1 + sum(arms[:2])
