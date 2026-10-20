"""Heteroacyclic parent hydrides (silane, disilane, disiloxane, phosphane, ...)
carrying -COOH, -CONH2, -CN or -CHO directly on a skeletal atom, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-65.1.2.2.2 / P-66.1.1.1.1.3 / P-66.5.1.1.3 / P-66.6.1.1.3: a group on a
  heteroacyclic parent hydride is always expressed by the suffix
  'carboxylic acid', 'carboxamide', 'carbonitrile' or 'carbaldehyde'
  ('disiloxanecarboxylic acid', 'phosphanecarboxamide').
- P-44.1.1: the parent hydride holding the principal characteristic group is
  the senior parent, so every other substituent becomes a prefix.
- Numbering: the suffix takes the lowest locants, then the prefixes
  (P-31.1.4.3.4); P-14.3.4.3 omits the locant of a sole substituent on a
  symmetrical parent, and P-14.3.4.2(a) on a mononuclear one.
"""

import contextlib

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    multiplied_word,
    nonstandard_bonding,
    specified_stereo_elements,
    substituent_locant_set_and_citation,
)
from ._hetero_prefixes import MONONUCLEAR_HYDRIDES
from ._hydride_chain import _STEMS, _chain_atoms
from ._numerals import multiplying_prefix
from ._substituents import cited_branch_stereo, format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_NITROGEN = 7
_CLASSES = ("acid", "amide", "amidine", "nitrile", "aldehyde")
_SUFFIX = {
    "acid": "carboxylic acid",
    "amide": "carboxamide",
    "amidine": "carboximidamide",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
}


def _group_class(mol, carbon):
    """The class of a terminal -COOH/-CONH2/-CN/-CHO carbon, else None."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.IsInRing() or atom.GetFormalCharge():
        return None
    neighbors = list(atom.GetNeighbors())
    order = {n.GetIdx(): mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() for n in neighbors}
    if any(order[n.GetIdx()] == 3.0 and n.GetAtomicNum() == _NITROGEN for n in neighbors):
        return "nitrile" if len(neighbors) == 2 else None
    imino = [
        n
        for n in neighbors
        if order[n.GetIdx()] == 2.0 and n.GetAtomicNum() == _NITROGEN and n.GetDegree() == 1 and n.GetTotalNumHs() == 1
    ]
    if len(imino) == 1:
        amino = [
            n
            for n in neighbors
            if n.GetIdx() != imino[0].GetIdx() and order[n.GetIdx()] == 1.0 and n.GetAtomicNum() == _NITROGEN
            and n.GetDegree() == 1 and n.GetTotalNumHs() == 2
        ]
        return "amidine" if len(amino) == 1 and len(neighbors) == 3 else None
    doubled = [n for n in neighbors if order[n.GetIdx()] == 2.0 and n.GetAtomicNum() == 8 and n.GetDegree() == 1]
    if len(doubled) != 1:
        return None
    others = [n for n in neighbors if n.GetIdx() != doubled[0].GetIdx()]
    if len(others) == 1 and atom.GetTotalNumHs() == 1:
        return "aldehyde"
    if len(others) != 2:
        return None
    hetero = [
        n
        for n in others
        if order[n.GetIdx()] == 1.0
        and (
            (n.GetDegree() == 1 and (n.GetAtomicNum(), n.GetTotalNumHs()) in ((8, 1), (_NITROGEN, 2)))
            or _substituted_amide_nitrogen(n, atom)
        )
    ]
    if len(hetero) != 1:
        return None
    return "acid" if hetero[0].GetAtomicNum() == 8 else "amide"


def _substituted_amide_nitrogen(nitrogen, carbon):
    """An amino nitrogen of the amide group whose substituents are carbon atoms that are not acyl, imidoyl or
    thioacyl carbons (P-66.1.1.3.1)."""
    if nitrogen.GetAtomicNum() != _NITROGEN or nitrogen.GetFormalCharge() or nitrogen.IsInRing() or nitrogen.GetDegree() < 2:
        return False
    return all(
        n.GetAtomicNum() == 6
        and not any(b.GetBondTypeAsDouble() != 1.0 and b.GetOtherAtom(n).GetAtomicNum() in (7, 8, 16, 34, 52) for b in n.GetBonds())
        for n in nitrogen.GetNeighbors()
        if n.GetIdx() != carbon.GetIdx()
    ) and all(b.GetBondTypeAsDouble() == 1.0 for b in nitrogen.GetBonds())


def _siloxane_skeleton(mol, graph):
    """The Si-O-Si-...-Si backbone (only bridging oxygens), or None."""
    bridges = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 8 and a.GetDegree() == 2 and all(n.GetAtomicNum() == 14 for n in a.GetNeighbors())
    ]
    if not bridges:
        return None
    silicon = {n for o in bridges for n in graph[o]}
    if len(silicon) != len(bridges) + 1 or any(sum(n in bridges for n in graph[s]) > 2 for s in silicon):
        return None
    ends = [s for s in silicon if sum(n in bridges for n in graph[s]) == 1]
    if len(ends) != 2:
        return None
    chain, previous = [ends[0]], None
    while True:
        step = [n for n in graph[chain[-1]] if n in bridges and n != previous]
        if not step:
            break
        bridge = step[0]
        (nxt,) = [n for n in graph[bridge] if n != chain[-1]]
        chain += [bridge, nxt]
        previous = bridge
    return chain if len(chain) == len(silicon) + len(bridges) else None


def _group_nitrogens(mol):
    return frozenset(
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == _NITROGEN and a.GetDegree() == 1 and all(n.GetAtomicNum() == 6 for n in a.GetNeighbors())
    )


def _separate_chains(mol, graph):
    """Parent hydrides for a molecule whose heteroatoms form several chains: every chain that is a parent hydride on its
    own, the others being substituents of the one that carries the principal group (P-44.1.1)."""
    ignore = _group_nitrogens(mol)
    pool = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in _STEMS and a.GetIdx() not in ignore and not a.IsInRing() and not a.GetFormalCharge()
    }
    parents, seen = [], set()
    for start in sorted(pool):
        if start in seen:
            continue
        component, stack = {start}, [start]
        while stack:
            for n in graph[stack.pop()]:
                if n in pool and n not in component:
                    component.add(n)
                    stack.append(n)
        seen |= component
        found = _chain_atoms(mol, graph, allow_double=True, ignore=frozenset(range(mol.GetNumAtoms())) - component)
        if found is not None:
            z, chain = found
            parents.append(("chain", chain, f"{multiplying_prefix(len(chain))}{_STEMS[z]}"))
    return parents


def _parent(mol, graph):
    """(kind, skeleton atoms, parent hydride name) of the heteroacyclic parent, or None."""
    found = _chain_atoms(mol, graph, skip_nitrogen=True) or _chain_atoms(mol, graph, allow_double=True)
    if found is None:
        group_nitrogens = _group_nitrogens(mol)
        found = _chain_atoms(mol, graph, allow_double=True, ignore=group_nitrogens) if group_nitrogens else None
    if found is not None:
        z, chain = found
        return "chain", chain, f"{multiplying_prefix(len(chain))}{_STEMS[z]}"
    siloxane = _siloxane_skeleton(mol, graph)
    if siloxane is not None:
        return "chain", siloxane, f"{multiplying_prefix((len(siloxane) + 1) // 2)}siloxane"
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in MONONUCLEAR_HYDRIDES]
    if len(centers) == 1 and not centers[0].GetFormalCharge() and not centers[0].IsInRing():
        return "mononuclear", [centers[0].GetIdx()], MONONUCLEAR_HYDRIDES[centers[0].GetAtomicNum()][0]
    return None


def _match(mol):
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    parent = _parent(mol, graph)
    candidates = [parent] if parent is not None else _separate_chains(mol, graph)
    matches = [m for m in (_match_parent(mol, graph, c) for c in candidates) if m is not None]
    return matches[0] if len(matches) == 1 else None


def _match_parent(mol, graph, parent):
    _, skeleton, _ = parent
    skeleton_set = set(skeleton)
    attached = {}
    for atom in skeleton:
        for n in graph[atom]:
            if n not in skeleton_set and mol.GetAtomWithIdx(n).GetAtomicNum() == 6:
                found = _group_class(mol, n)
                if found is not None:
                    attached[n] = (atom, found)
    if not attached:
        return None
    principal = next(c for c in _CLASSES if any(cls == c for _, cls in attached.values()))
    attached_principal = {c: host for c, (host, cls) in attached.items() if cls == principal}
    elsewhere = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in attached and _group_class(mol, a.GetIdx()) is not None
    ]
    senior = [i for i in elsewhere if _CLASSES.index(_group_class(mol, i)) <= _CLASSES.index(principal)]
    if senior and not _outnumbered_by_skeleton(mol, skeleton_set, senior, len(attached_principal)):
        return None
    return parent, attached_principal, principal


def _outnumbered_by_skeleton(mol, skeleton_set, senior, count):
    """Whether no substituent branch holds more groups of the principal class than the skeleton: with equal numbers the
    heteroatom parent hydride is senior to the carbon parent (P-44.1.1, P-44.1.2.1)."""
    branch_of = {}
    for start in senior:
        if start in branch_of:
            continue
        seen, stack = {start}, [start]
        while stack:
            for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
                if n.GetIdx() not in skeleton_set and n.GetIdx() not in seen:
                    seen.add(n.GetIdx())
                    stack.append(n.GetIdx())
        for a in seen:
            branch_of[a] = seen
    return all(sum(1 for i in senior if i in branch_of[c]) <= count for c in senior)


def _is_oxoacid_derivative(mol) -> bool:
    """A Group 15 atom with a double-bonded chalcogen is an oxoacid centre: its cyanides are named as acid derivatives."""
    return any(
        a.GetAtomicNum() in (15, 33, 51)
        and any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(a).GetAtomicNum() in (8, 16) for b in a.GetBonds()
        )
        for a in mol.GetAtoms()
    )


def has_hydride_carbo_suffix_shape(mol) -> bool:
    return not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) and _match(mol) is not None and not _is_oxoacid_derivative(mol)


def name_hydride_carbo_suffix(mol) -> str:
    (kind, skeleton, parent_name), principal_carbons, principal = _match(mol)
    stereo = specified_stereo_elements(mol)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    skeleton_set = set(skeleton)
    word = _SUFFIX[principal]
    count = len(principal_carbons)
    bond_codes = _skeleton_bond_descriptors(mol, stereo, skeleton_set, set(principal_carbons))

    entries = []
    nitrogen_entries = {}
    nitrogens = {c: _amide_nitrogen(mol, c, host) for c, host in principal_carbons.items()} if principal == "amide" else {}
    group_atoms = set(principal_carbons) | {n for n in nitrogens.values() if n is not None}
    roots = [n for atom in skeleton for n in graph[atom] if n not in skeleton_set and n not in principal_carbons]
    roots += [m for n in nitrogens.values() if n is not None for m in graph[n] if m not in group_atoms]
    with cited_branch_stereo(mol, graph, skeleton_set | group_atoms, roots) if stereo else contextlib.nullcontext():
        for atom in skeleton:
            for n in graph[atom]:
                if n in skeleton_set or n in principal_carbons:
                    continue
                entries.append((atom, *name_branch(graph, n, atom, halogens, mol=mol, unsaturated=True)))
        for carbon, nitrogen in nitrogens.items():
            nitrogen_entries[carbon] = [
                name_branch(graph, m, nitrogen, halogens, mol=mol, unsaturated=True) for m in graph[nitrogen] if m != carbon
            ]

    if kind == "mononuclear":
        parent_name = _bonding_marked(_lambda_marks(mol, {skeleton[0]: 1}, 1), parent_name)
        if any(nitrogen_entries.values()):
            grouped = group_substituents(
                _merge(entries, {skeleton[0]: 1}) | _nitrogen_locants(nitrogen_entries, {c: 1 for c in nitrogen_entries})
            )
            return _join_parent(format_substituent_prefixes(grouped, omit_locants=True), parent_name) + multiplied_word(count, word)
        prefix = format_mononuclear_prefixes([(name, compound) for _, name, compound in entries]) if entries else ""
        return _join_parent(prefix, parent_name) + multiplied_word(count, word)

    best = None
    for candidate in (skeleton, skeleton[::-1]):
        position_of = {a: i + 1 for i, a in enumerate(candidate)}
        suffix_locants = sorted(position_of[host] for host in principal_carbons.values())
        ene_locants = [
            i + 1 for i, (a, b) in enumerate(zip(candidate, candidate[1:])) if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() == 2.0
        ]
        grouped = group_substituents(_merge(entries, position_of))
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        on_nitrogen = _nitrogen_locants(nitrogen_entries, {c: position_of[host] for c, host in principal_carbons.items()})
        descriptors = [
            (i + 1, bond_codes[frozenset((a, b))])
            for i, (a, b) in enumerate(zip(candidate, candidate[1:]))
            if frozenset((a, b)) in bond_codes
        ]
        marks = _lambda_marks(mol, position_of, len(skeleton))
        key = (suffix_locants, ene_locants, locant_set, citation, [locant for locant, code in descriptors if code == "Z"], marks)
        if best is None or key < best[0]:
            best = (key, suffix_locants, ene_locants, group_substituents(_merge(entries, position_of) | on_nitrogen), descriptors, marks)

    _, suffix_locants, ene_locants, grouped, descriptors, marks = best
    total = count + sum(len(info["locants"]) for info in grouped.values())
    siloxane = parent_name.endswith("siloxane")
    diazene_like = len(skeleton) == 2 and bool(ene_locants)
    symmetric_sole = (total == 1 and len(skeleton) == (3 if siloxane else 2) and not ene_locants) or diazene_like
    body = multiplied_word(count, word)
    prefix = format_substituent_prefixes(grouped, omit_all=diazene_like and not any(nitrogen_entries.values()))
    if ene_locants:
        stem = parent_name[:-3]
        if len(skeleton) == 2:
            parent_name = f"{stem}ene"
        else:
            parent_name = f"{stem}{'a' if len(ene_locants) > 1 else ''}-{','.join(map(str, ene_locants))}-{multiplied_word(len(ene_locants), 'ene')}"
    cited = ",".join(code if len(skeleton) == 2 else f"{locant}{code}" for locant, code in descriptors)
    lead = f"({cited})-" if cited else ""
    parent_text = _join_parent(prefix, _bonding_marked(marks, parent_name))
    if symmetric_sole:
        return f"{lead}{parent_text}{body}"
    return f"{lead}{parent_text}-{','.join(map(str, suffix_locants))}-{body}"


def _lambda_marks(mol, position_of, size):
    """'1λ5,2λ5' (or 'λ5' for a mononuclear parent) for the skeletal atoms of nonstandard bonding number (P-14.1)."""
    marks = [
        (position, nonstandard_bonding(mol.GetAtomWithIdx(atom)))
        for atom, position in sorted(position_of.items(), key=lambda item: item[1])
    ]
    return tuple((position if size > 1 else None, bonding) for position, bonding in marks if bonding)


def _bonding_marked(marks, parent_name):
    """`parent_name` led by its lambda marks, e.g. 'λ5-phosphane' (None when no atom is nonstandard)."""
    if not marks:
        return (None, parent_name)
    return (",".join(f"{'' if position is None else position}λ{bonding}" for position, bonding in marks), parent_name)


def _join_parent(prefix, marked):
    """Prefixes, then the lambda marks hyphenated from both neighbours, then the parent hydride name."""
    marks, parent_name = marked
    if marks is None:
        return prefix + parent_name
    return f"{prefix}{'-' if prefix else ''}{marks}-{parent_name}"


def _amide_nitrogen(mol, carbon, host):
    """The nitrogen of the amide group on `carbon`, or None when the group is not an amide."""
    return next(
        (n.GetIdx() for n in mol.GetAtomWithIdx(carbon).GetNeighbors() if n.GetAtomicNum() == _NITROGEN and n.GetIdx() != host),
        None,
    )


def _nitrogen_locants(nitrogen_entries, host_position):
    """{N, N', ...: substituents} of the amide nitrogens, the primes following the order of the suffix locants
    (P-66.1.1.3.1)."""
    ordered = sorted(nitrogen_entries, key=lambda carbon: host_position[carbon])
    return {"N" + "'" * i: nitrogen_entries[carbon] for i, carbon in enumerate(ordered) if nitrogen_entries[carbon]}


def _skeleton_bond_descriptors(mol, stereo, skeleton_set, group_carbons):
    """{bond atom pair: descriptor} of the stereogenic skeleton double bonds, cited with the parent (P-93.4.2.1.1); a
    stereo element of the parent atoms or the group carbons that is not such a bond is not supported."""
    codes = {}
    parent_atoms = skeleton_set | group_carbons
    for kind, index, code in stereo or []:
        if kind == "bond":
            bond = mol.GetBondWithIdx(index)
            ends = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            if all(end in skeleton_set for end in ends):
                codes[frozenset(ends)] = code
            elif any(end in parent_atoms for end in ends):
                raise UnsupportedStructure("this stereo element of a heteroacyclic parent hydride is not supported yet")
        elif index in parent_atoms:
            raise UnsupportedStructure("this stereo element of a heteroacyclic parent hydride is not supported yet")
    return codes


def _merge(entries, position_of):
    merged = {}
    for atom, name, compound in entries:
        merged.setdefault(position_of[atom], []).append((name, compound))
    return merged
