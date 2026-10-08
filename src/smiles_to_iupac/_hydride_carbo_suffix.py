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

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    multiplied_word,
    specified_stereo_elements,
    substituent_locant_set_and_citation,
)
from ._hetero_prefixes import MONONUCLEAR_HYDRIDES
from ._hydride_chain import _STEMS, _chain_atoms
from ._numerals import multiplying_prefix
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

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
        if n.GetDegree() == 1
        and order[n.GetIdx()] == 1.0
        and (n.GetAtomicNum(), n.GetTotalNumHs()) in ((8, 1), (_NITROGEN, 2))
    ]
    if len(hetero) != 1:
        return None
    return "acid" if hetero[0].GetAtomicNum() == 8 else "amide"


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


def _parent(mol, graph):
    """(kind, skeleton atoms, parent hydride name) of the heteroacyclic parent, or None."""
    if mol.GetRingInfo().NumRings():
        return None
    found = _chain_atoms(mol, graph, skip_nitrogen=True)
    if found is not None:
        z, chain = found
        return "chain", chain, f"{multiplying_prefix(len(chain))}{_STEMS[z]}"
    siloxane = _siloxane_skeleton(mol, graph)
    if siloxane is not None:
        return "chain", siloxane, f"{multiplying_prefix((len(siloxane) + 1) // 2)}siloxane"
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in MONONUCLEAR_HYDRIDES]
    if len(centers) == 1 and not centers[0].GetFormalCharge():
        return "mononuclear", [centers[0].GetIdx()], MONONUCLEAR_HYDRIDES[centers[0].GetAtomicNum()][0]
    return None


def _match(mol):
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    parent = _parent(mol, graph)
    if parent is None:
        return None
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
    elsewhere = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in attached and _group_class(mol, a.GetIdx()) is not None
    ]
    if any(_CLASSES.index(_group_class(mol, i)) <= _CLASSES.index(principal) for i in elsewhere):
        return None
    return parent, {c: host for c, (host, cls) in attached.items() if cls == principal}, principal


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
    return _match(mol) is not None and not _is_oxoacid_derivative(mol)


def name_hydride_carbo_suffix(mol) -> str:
    (kind, skeleton, parent_name), principal_carbons, principal = _match(mol)
    if specified_stereo_elements(mol):
        raise UnsupportedStructure("stereodescriptors on a heteroacyclic parent hydride are not supported yet")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    skeleton_set = set(skeleton)
    word = _SUFFIX[principal]
    count = len(principal_carbons)

    entries = []
    for atom in skeleton:
        for n in graph[atom]:
            if n in skeleton_set or n in principal_carbons:
                continue
            entries.append((atom, *name_branch(graph, n, atom, halogens, mol=mol, unsaturated=True)))

    if kind == "mononuclear":
        prefix = format_mononuclear_prefixes([(name, compound) for _, name, compound in entries]) if entries else ""
        return f"{prefix}{parent_name}{multiplied_word(count, word)}"

    best = None
    for candidate in (skeleton, skeleton[::-1]):
        position_of = {a: i + 1 for i, a in enumerate(candidate)}
        suffix_locants = sorted(position_of[host] for host in principal_carbons.values())
        grouped = group_substituents(_merge(entries, position_of))
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        key = (suffix_locants, locant_set, citation)
        if best is None or key < best[0]:
            best = (key, suffix_locants, grouped)

    _, suffix_locants, grouped = best
    total = count + sum(len(info["locants"]) for info in grouped.values())
    siloxane = parent_name.endswith("siloxane")
    symmetric_sole = total == 1 and len(skeleton) == (3 if siloxane else 2)
    body = multiplied_word(count, word)
    prefix = format_substituent_prefixes(grouped)
    if symmetric_sole:
        return f"{prefix}{parent_name}{body}"
    return f"{prefix}{parent_name}-{','.join(map(str, suffix_locants))}-{body}"


def _merge(entries, position_of):
    merged = {}
    for atom, name, compound in entries:
        merged.setdefault(position_of[atom], []).append((name, compound))
    return merged
