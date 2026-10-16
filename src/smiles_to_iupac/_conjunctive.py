"""Conjunctive names (P-13.5.2, P-15.6): the cyclic component and the side chain that carries the principal characteristic
group are juxtaposed, as in benzeneacetic acid and naphthalene-2-propanol. The ring is numbered with its points of
attachment as principal positions, the side chain is named as the unbranched parent compound, and the substituents of
both are cited in front, those of the side chain under Greek letters counted from the principal group."""

import re
from collections import namedtuple

from rdkit import Chem

from ._common import (
    KEEP_ASSEMBLY_ENDING, RuleUndefinedStructure, UnsupportedStructure, adjacency, alpha_sort_key, assembly_join, assign_cip_labels,
    group_substituents, halogen_substituents,
)
from ._diester_ring_diyl import PARENT_START, _stereo_locant, _system_of, evaluate_skeleton
from ._numerals import multiplying_prefix
from ._parent_locants import locant_key
from ._polyfunctional import assembly_substituent
from ._ring_system_seniority import ring_seniority_key
from ._substituents import cited_branch_stereo, format_substituent_prefixes, name_branch

_GREEK = "αβγδεζηθικλμνξοπρστυφχψω"
_Class = namedtuple("_Class", "name query offset probe")
_CLASSES = (
    _Class("acid", Chem.MolFromSmarts("[CX3](=O)[OX2H1]"), 1, lambda n: "C" * (n - 1) + "C(=O)O"),
    _Class("amide", Chem.MolFromSmarts("[CX3](=O)[NX3H2]"), 1, lambda n: "C" * (n - 1) + "C(N)=O"),
    _Class("nitrile", Chem.MolFromSmarts("[CX2]#[NX1]"), 1, lambda n: "C" * (n - 1) + "C#N"),
    _Class("aldehyde", Chem.MolFromSmarts("[CX3H1](=O)[#6]"), 1, lambda n: "C" * (n - 1) + "C=O"),
    _Class("ketone", Chem.MolFromSmarts("[#6][CX3](=O)[#6]"), 0, None),
    _Class("alcohol", Chem.MolFromSmarts("[#6][OX2H1]"), 0, lambda n: "C" * n + "O"),
    _Class("thiol", Chem.MolFromSmarts("[#6][SX2H1]"), 0, lambda n: "C" * n + "S"),
    _Class("amine", Chem.MolFromSmarts("[#6][NX3H2]"), 0, lambda n: "C" * n + "N"),
)
_DIRECT = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in {
        "acid": "[R][CX3](=O)[OX2H1]", "amide": "[R][CX3](=O)[NX3H2]", "nitrile": "[R][CX2]#[NX1]",
        "aldehyde": "[R][CX3H1]=O", "alcohol": "[R;#6][OX2H1]", "thiol": "[R;#6][SX2H1]", "amine": "[R;#6][NX3H2]",
    }.items()
}
_SENIOR_OTHERS = tuple(
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=[O,S,Se,Te])[O,S,Se,Te,N,F,Cl,Br,I]",
        "[CX3](=[NX2])[O,S,Se,Te,N]",
        "[#6X2](=*)=*",
        "[#15,#16,#34,#52,#33,#51,#5,#14;$(*~[OX1])]-[#8,#7,#9,#17,#35,#53]",
    )
)
_OTHER_AMINE = Chem.MolFromSmarts("[NX3;!R;!H2;!$(N-[!#6;!#1]);!$(N-C=[O,S,N]);!$(N-a)]")
_NITRO = Chem.MolFromSmarts("[NX3+](=O)[O-]")
_Chain = namedtuple("_Chain", "path ring_atom rings atoms")
_MARK = "\u00a7"
_SUPERSCRIPT = "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079"
_ORGANIC = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 52, 53}
_ASSEMBLY_TAIL = re.compile(r"(?P<core>.+)-(?P<spot>\d+[a-z]*\u2032*(?:\([^)]*\))?)-yl")
_DESCRIPTOR_LEAD = re.compile(r"\([^()]*\)-")
_TAIL = re.compile(r"(?P<core>.+?)-(?P<locants>[\d,]+(?:\([^)]*\))?)-(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?\u00a7")


class _Mark(str):
    """A Greek letter locant, with the superscript number of the side chain when several are identical, ordered by letter
    and then by number."""

    def __new__(cls, position, number):
        mark = super().__new__(cls, _GREEK[position] + (_SUPERSCRIPT[number] if number else ""))
        mark.order = (position, number)
        return mark

    def __lt__(self, other):
        return self.order < other.order

    def __gt__(self, other):
        return self.order > other.order


def conjunctive_names(mol, name_parent):
    """[(name, rule)] with the conjunctive name of the structure, or [] when it has no side chain with a principal group
    on a ring; raises UnsupportedStructure when it has one this module cannot name."""
    nitro = {i for match in mol.GetSubstructMatches(_NITRO) for i in match}
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        (a.GetFormalCharge() and a.GetIdx() not in nitro)
        or a.GetNumRadicalElectrons()
        or a.GetIsotope()
        or a.GetAtomicNum() not in _ORGANIC
        for a in mol.GetAtoms()
    ):
        return []
    klass = _principal_class(mol) if mol.GetRingInfo().NumRings() else None
    if klass is None:
        return []
    graph = adjacency(mol)
    anchors = {m[0] for m in mol.GetSubstructMatches(klass.query)}
    chains = [chain for anchor in sorted(anchors) if (chain := _side_chain(mol, graph, klass, anchor)) is not None]
    if not chains:
        return []
    options = _chosen_chains(mol, graph, klass, anchors, chains)
    if options is None:
        return []
    built = [_assemble(mol, graph, klass, anchors, selected, name_parent) for selected in options]
    best = min(built, key=lambda b: b[0])
    if len({name for key, name in built if key == best[0]}) > 1:
        raise UnsupportedStructure("several side chains are equally senior as the conjunctive parent")
    return [(best[1], "P-15.6.1.1")]


def _principal_class(mol):
    """The class of the principal characteristic group when a conjunctive name can be built for it, else None."""
    plain = {frozenset(m) for klass in _CLASSES[:2] for m in mol.GetSubstructMatches(klass.query)}
    if any(frozenset(m) not in plain for pattern in _SENIOR_OTHERS for m in mol.GetSubstructMatches(pattern)):
        return None
    klass = next((k for k in _CLASSES if mol.HasSubstructMatch(k.query)), None)
    if klass is None or klass.probe is None:
        return None
    if klass.name == "amine" and mol.HasSubstructMatch(_OTHER_AMINE):
        return None
    return klass


def _side_chain(mol, graph, klass, anchor):
    """The side chain from the principal group at `anchor` to the ring it ends at, running to the senior ring or ring
    system (P-15.6.1.1, P-15.6.1.3, P-15.6.1.6)."""
    if mol.GetAtomWithIdx(anchor).IsInRing():
        return None
    acyclic = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
    parent, stack = {anchor: None}, [anchor]
    while stack:
        current = stack.pop()
        for n in graph[current]:
            if n in acyclic and n not in parent:
                parent[n] = current
                stack.append(n)
    candidates = []
    for x in parent:
        for r in graph[x]:
            if not mol.GetAtomWithIdx(r).IsInRing():
                continue
            path = [x]
            while parent[path[-1]] is not None:
                path.append(parent[path[-1]])
            path.reverse()
            bond = mol.GetBondBetweenAtoms(x, r)
            if not _saturated(mol, path, klass) or (len(path) == 1 and klass.offset):
                continue
            if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
                continue
            rings, atoms = _system_of(mol, r)
            candidates.append((ring_seniority_key(mol, atoms), -len(path), x, r, path, rings, atoms))
    if not candidates:
        return None
    candidates.sort(key=lambda c: c[:2])
    best = [c for c in candidates if c[:2] == candidates[0][:2]]
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    if len({(ranks[c[2]], ranks[c[3]]) for c in best}) > 1:
        raise UnsupportedStructure("several rings are equally senior as the cyclic component")
    _, _, _, ring_atom, path, rings, atoms = best[0]
    return _Chain(path, ring_atom, rings, atoms)


def _saturated(mol, path, klass):
    """Only carbon atoms joined by single bonds, the principal group aside (P-15.6.2 (a), (b))."""
    for i, index in enumerate(path):
        atom = mol.GetAtomWithIdx(index)
        if i == 0 and klass.offset:
            continue
        if atom.GetAtomicNum() != 6 or atom.GetHybridization() != Chem.HybridizationType.SP3:
            return False
    return all(mol.GetBondBetweenAtoms(a, b).GetBondType() == Chem.BondType.SINGLE for a, b in zip(path, path[1:]))


def _reach(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n not in blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _chain_substituents(mol, graph, klass, chain, anchors):
    """[(Greek position, name, is_compound)] of what the side chain carries (P-15.6.1.2, P-15.6.1.3), or None when another
    principal group among them rules the chain out (P-15.6.2 (c))."""
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    path = chain.path
    owned = {i for m in mol.GetSubstructMatches(klass.query) if m[0] == path[0] for i in m[1:]}
    found = []
    for i, index in enumerate(path):
        for n in graph[index]:
            if n in path or n in owned or (index == path[-1] and n == chain.ring_atom):
                continue
            if anchors & _reach(graph, n, set(path)):
                return None
            if i - klass.offset >= len(_GREEK):
                raise RuleUndefinedStructure("the side chain is longer than the Greek alphabet (P-15.6.1.2)")
            name, compound = name_branch(graph, n, index, halogens, aromatic, mol=mol)
            found.append((i - klass.offset, name, compound))
    return found


def _chosen_chains(mol, graph, klass, anchors, chains):
    """The identical side chains of one cyclic component that express the most principal groups, then the senior ring and
    the longer chain (P-15.6.1.5, P-15.6.1.6); several selections are returned when these criteria leave a tie."""
    groups = {}
    for chain in chains:
        found = _chain_substituents(mol, graph, klass, chain, anchors)
        if found is not None:
            subs = tuple(sorted(found))
            groups.setdefault((frozenset(chain.atoms), len(chain.path), subs), []).append((chain, subs))
    if not groups:
        return None

    def rank(item):
        return -len(item[1]), ring_seniority_key(mol, item[0][0]), -item[0][1]

    ranked = sorted(groups.items(), key=rank)
    return [group[1] for group in ranked if rank(group) == rank(ranked[0])]


def _assembly_partner(mol, graph, atoms):
    """The atoms of the identical ring system joined to `atoms` by a single bond, which together form a ring assembly
    (P-28), else None."""
    skeleton = Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(atoms), canonical=True)
    for a in atoms:
        for n in graph[a]:
            if n not in atoms and mol.GetAtomWithIdx(n).IsInRing():
                other = _system_of(mol, n)[1]
                if Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(other), canonical=True) == skeleton:
                    return other
    return None


def _locant_text(locant):
    return re.sub(r"^(\d+)", lambda m: m.group(1) + "\u2032" * locant.primes, str(locant.number))


def _assemble(mol, graph, klass, anchors, selected, name_parent):
    """(key, name): the key ranks alternative selections by the prefixes of the whole name (P-45.2.1 to P-45.4)."""
    chains = [chain for chain, _ in selected]
    skeleton = set(chains[0].atoms) | {a for chain in chains for a in chain.path}
    roots = [n for a in skeleton for n in graph[a] if n not in skeleton]
    with cited_branch_stereo(mol, graph, skeleton, roots):
        return _assembled(mol, graph, klass, anchors, chains, name_parent)


def _assembled(mol, graph, klass, anchors, chains, name_parent):
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    rings, atoms = chains[0].rings, chains[0].atoms
    blocked = set()
    for chain in chains:
        branch = _reach(graph, chain.path[-1], atoms)
        if any(n in atoms and (a, n) != (chain.path[-1], chain.ring_atom) for a in branch for n in graph[a]):
            raise UnsupportedStructure("a side chain is joined to the ring system more than once")
        blocked |= branch
    count = len(chains)
    attach = [chain.ring_atom for chain in chains]
    if sum(m[0] in atoms for m in mol.GetSubstructMatches(_DIRECT[klass.name])) > 1:
        raise UnsupportedStructure("P-15.6.2 (e): several groups of the principal kind on the ring system")
    partner = _assembly_partner(mol, graph, atoms)
    if partner is None:
        grouped, core, locants, position_of, ring_labels, saturated_ring = _ring_part(
            mol, graph, rings, atoms, attach, blocked, halogens, aromatic
        )
    else:
        if count > 1 or any(
            (other := _side_chain(mol, graph, klass, anchor)) is not None and other.ring_atom in partner for anchor in anchors
        ):
            raise UnsupportedStructure("several side chains on a ring assembly")
        grouped, core, locants, position_of, ring_labels = _assembly_part(
            mol, graph, chains[0], atoms | partner, blocked, halogens, aromatic
        )
        saturated_ring = False
    order = sorted(range(count), key=lambda k: position_of[chains[k].ring_atom])
    if count > len(_SUPERSCRIPT) - 1:
        raise UnsupportedStructure("too many identical side chains for superscript locants")
    probe = Chem.Mol(mol)
    assign_cip_labels(probe)
    chain_labels = []
    for rank, k in enumerate(order, 1):
        for position, name, compound in sorted(_chain_substituents(mol, graph, klass, chains[k], anchors)):
            mark = _Mark(position, rank if count > 1 else 0)
            grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(mark)
        for i, index in enumerate(chains[k].path):
            if i >= klass.offset and probe.GetAtomWithIdx(index).HasProp("_CIPCode"):
                chain_labels.append((_Mark(i - klass.offset, rank if count > 1 else 0), probe.GetAtomWithIdx(index).GetProp("_CIPCode")))
    _require_every_stereo_element_cited(mol, atoms | {a for chain in chains for a in chain.path}, blocked, graph)
    descriptor = ""
    if chain_labels or ring_labels:
        cited = [f"{mark}{code}" for mark, code in sorted(chain_labels, key=lambda item: item[0].order)]
        cited += [f"{locant}{code}" for locant, code in ring_labels]
        descriptor = f"({','.join(cited)})-"
    prefixes = format_substituent_prefixes(grouped)
    cited = "" if (core == "benzene" or saturated_ring) and count == 1 else f"-{locants}-"
    chain = re.sub(r"-1-(?=ol$|thiol$|amine$)", "", name_parent(klass.probe(len(chains[0].path))))
    side = chain if count == 1 else _multiplied(count, chain, len(chains[0].path))
    if partner is None:
        body = f"{prefixes}{'-' if prefixes and core[0].isdigit() else ''}{core}{cited}{side}"
    else:
        body = assembly_join(prefixes, f"{core}{cited}{side}")
    return _prefix_key(grouped), descriptor + body


def _ring_part(mol, graph, rings, atoms, attach, blocked, halogens, aromatic):
    """(prefix groups, cyclic component name, attachment locants, positions, stereo labels, saturated monocycle)."""
    found = evaluate_skeleton(mol, graph, "ring", rings, atoms, attach, blocked, _MARK)
    if found is None:
        raise UnsupportedStructure("this cyclic component has no supported numbering")
    group_name, position_of, ring_stereo = found[1], found[2], found[3]
    ring_labels = []
    if ring_stereo:
        if any(code.startswith("(") for _, code in ring_stereo):
            raise UnsupportedStructure("this stereodescriptor of the cyclic component is not cited with a locant")
        group_name = group_name[_DESCRIPTOR_LEAD.match(group_name).end():]
        ring_labels = sorted(
            ((_stereo_locant(a, position_of), code) for a, code in ring_stereo), key=lambda item: locant_key(item[0])
        )
    start = PARENT_START.get()
    ring_subs = {}
    for a in atoms:
        for n in graph[a]:
            if n not in atoms and n not in blocked:
                ring_subs.setdefault(position_of[a], []).append(name_branch(graph, n, a, halogens, aromatic, mol=mol))
    grouped = group_substituents(ring_subs)
    if format_substituent_prefixes(grouped) != group_name[:start]:
        raise UnsupportedStructure("the prefixes of the cyclic component are not those of its group name")
    benzene = len(rings) == 1 and len(atoms) == 6 and all(
        mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in atoms
    )
    tail = _TAIL.fullmatch(group_name[start:].lstrip("-"))
    if tail is None and not benzene:
        raise UnsupportedStructure("the group name of the cyclic component has no recognized attachment locants")
    core = "benzene" if benzene else tail.group("core")
    locants = ",".join(str(x) for x in sorted(position_of[a] for a in attach)) if benzene else tail.group("locants")
    saturated_ring = len(rings) == 1 and all(
        mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and mol.GetAtomWithIdx(a).GetHybridization() == Chem.HybridizationType.SP3
        for a in atoms
    )
    return grouped, core, locants, position_of, ring_labels, saturated_ring


def _assembly_part(mol, graph, chain, assembly_atoms, blocked, halogens, aromatic):
    """The same parts for a ring assembly, read from the name of the assembly as a substituent group entered at the ring
    atom the side chain ends at (P-28, P-15.6.1.2)."""
    for element in Chem.FindPotentialStereo(mol):
        if element.specified == Chem.StereoSpecified.Specified and element.centeredOn in assembly_atoms:
            raise UnsupportedStructure("a stereo element of a ring assembly is not cited with a locant")
    place = {}
    token = KEEP_ASSEMBLY_ENDING.set(True)
    try:
        named = assembly_substituent(
            mol, graph, chain.ring_atom, chain.path[-1], halogens, aromatic, place=place, place_monocycles=True
        )
    finally:
        KEEP_ASSEMBLY_ENDING.reset(token)
    if named is None or "locants" not in place:
        raise UnsupportedStructure("this ring assembly has no supported group name")
    name = named[0].replace("'", "\u2032")
    core = name[place["start"]:]
    tail = _ASSEMBLY_TAIL.fullmatch(core)
    if tail is None:
        raise UnsupportedStructure("the group name of the ring assembly has no recognized attachment locant")
    position_of = {a: _locant_text(loc) for a, loc in place["locants"].items()}
    ring_subs = {}
    for a in assembly_atoms:
        for n in graph[a]:
            if n not in assembly_atoms and n not in blocked:
                ring_subs.setdefault(position_of[a], []).append(name_branch(graph, n, a, halogens, aromatic, mol=mol))
    grouped = group_substituents(ring_subs)
    if assembly_join(format_substituent_prefixes(grouped), core) != name:
        raise UnsupportedStructure("the prefixes of the ring assembly are not those of its group name")
    return grouped, tail.group("core"), tail.group("spot"), position_of, []


def _prefix_key(grouped):
    """P-45.2.1 more prefixes, P-45.2.2 lower locant set, P-45.2.3 lower locants in order of citation, P-45.4 earlier in
    alphanumerical order; a Greek letter locant of the side chain ranks before a numeral."""

    def order(locant):
        return (0, *locant.order) if isinstance(locant, _Mark) else (1, *locant_key(locant))

    names = sorted(grouped, key=alpha_sort_key)
    all_locants = sorted(order(x) for info in grouped.values() for x in info["locants"])
    cited = [sorted(order(x) for x in grouped[n]["locants"]) for n in names]
    return -len(all_locants), all_locants, cited, [alpha_sort_key(n) for n in names]


def _require_every_stereo_element_cited(mol, skeleton, blocked, graph):
    """A specified element must lie in the skeleton or wholly in a substituent, where a name cites it."""
    for element in Chem.FindPotentialStereo(mol):
        if element.specified != Chem.StereoSpecified.Specified:
            continue
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            continue
        bond = mol.GetBondWithIdx(element.centeredOn)
        ends = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
        if len(ends & skeleton) == 1:
            raise UnsupportedStructure("a stereodefined double bond joins the skeleton and a substituent")


def _multiplied(count, chain, length):
    """P-16.3 (c), (d): the name of a chain of five or more atoms begins with a multiplying prefix and is enclosed."""
    return f"{multiplying_prefix(count)}({chain})" if length >= 5 else multiplying_prefix(count) + chain
