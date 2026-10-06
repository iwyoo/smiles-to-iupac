"""Name assembly for a molecule matched on an Appendix 3 parent: unsaturation (P-101.6), principal group and prefixes
(P-101.7.1), esters and N-substituted amides/amines (P-65, P-66) and substituent groups (P-101.7.3)."""

from __future__ import annotations

import re
from dataclasses import dataclass

from rdkit import Chem

from ._appendix3_groups import SENIORITY
from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    multiplied_word,
)
from ._numerals import multiplying_prefix
from ._substituents import BRANCH_STEREO, format_substituent_prefixes, name_branch

SUFFIX = {
    ("acid", "c"): "carboxylic acid",
    ("acid", "o"): "oic acid",
    ("sulfonic", "c"): "sulfonic acid",
    ("anhydride", "c"): "carboxylic",
    ("anhydride", "o"): "oic",
    ("ester", "c"): "carboxylate",
    ("ester", "o"): "oate",
    ("halide", "c"): "carbonyl",
    ("halide", "o"): "oyl",
    ("amide", "c"): "carboxamide",
    ("amide", "o"): "amide",
    ("sulfonamide", "c"): "sulfonamide",
    ("amidine", "c"): "carboximidamide",
    ("amidine", "o"): "imidamide",
    ("hydrazide", "c"): "carbohydrazide",
    ("hydrazide", "o"): "hydrazide",
    ("nitrile", "c"): "carbonitrile",
    ("nitrile", "o"): "nitrile",
    ("aldehyde", "c"): "carbaldehyde",
    ("aldehyde", "o"): "al",
    ("ketone", ""): "one",
    ("alcohol", ""): "ol",
    ("thiol", ""): "thiol",
    ("amine", ""): "amine",
    ("imine", ""): "imine",
}
N_CLASSES = ("amide", "sulfonamide", "amidine", "hydrazide", "amine", "imine")
_INDICATED = re.compile(r"(\d+[a-c]?[¹²³⁴⁵⁶⁷⁸⁹]?′*)H-(.+)")
ENDINGS = {"anine": "enine", "ane": "ene", "an": "ene"}


@dataclass
class Choice:
    mapping: dict
    groups: list
    branches: list
    attach: tuple | None
    lost: list
    gained: list
    unsaturated_labels: frozenset
    hydro: list


def principal_class(groups):
    present = {g.cls for g in groups}
    for cls in SENIORITY:
        if cls == "ester":
            if "ester" in present:
                return "ester"
            if "ester_o" in present:
                return "ester_o"
        elif cls in present:
            return cls
    return None


def _elide(stem, word):
    return stem[:-1] if stem.endswith("e") and word[0] in "aeiouy" else stem


def _locant_text(label, faces, atom):
    return f"{label}{faces.get(atom, '')}"


def _ene_text(pair, sort_key):
    from ._appendix3_skeletons import _adjacent

    low, high = pair
    return low if _adjacent(pair) else f"{low}({high})"


def _words(stem, ene, yne, ending, principal_word, principal_locants):
    saturated = next(e for e in ENDINGS if ending == e)
    segments = []
    if ene:
        segments.append((ene, multiplied_word(len(ene), ENDINGS[saturated])))
    if yne:
        segments.append((yne, multiplied_word(len(yne), "yne")))
    if not segments:
        segments.append(([], saturated))
    if principal_word:
        segments.append((principal_locants, principal_word))
    text = stem + ("a" if segments[0][0] and segments[0][1].startswith(("di", "tri", "tetra", "penta", "hexa")) else "")
    for i, (locants, word) in enumerate(segments):
        if word.endswith("e") and i + 1 < len(segments) and segments[i + 1][1][0] in "aeiouy":
            word = word[:-1]
        text += (f"-{','.join(locants)}-" if locants else "") + word
    return text


def _without_added_hydrogen(hydro, skeleton, choice, members, principal, suffix_locants, sort_key):
    """P-31.1.4.3.4: a ketone suffix on an atom that the parent has doubly bonded takes an 'added hydrogen' on its
    partner, cited in brackets after the suffix locant; both leave the hydro prefixes. A parent that already has an
    indicated hydrogen moves it to the ketone atom when that atom has no partner (P-31.1.4.2.4, '2H-inden-2-one').
    Returns the remaining hydro atoms and the parent name."""
    remaining = list(hydro)
    name = skeleton.name
    indicated = _INDICATED.match(name)
    if principal == "ketone":
        neighbors = {}
        for bond in skeleton.query.GetBonds():
            a, b = skeleton.labels[bond.GetBeginAtomIdx()], skeleton.labels[bond.GetEndAtomIdx()]
            neighbors.setdefault(a, set()).add(b)
            neighbors.setdefault(b, set()).add(a)
        for index, group in enumerate(sorted(members, key=lambda g: sort_key(g.label))):
            if group.label not in remaining:
                continue
            partners = [p for p in neighbors.get(group.label, ()) if p in remaining and p != group.label]
            if not partners and indicated is not None and indicated.group(1) not in remaining:
                remaining.remove(group.label)
                remaining.append(indicated.group(1))
                name = f"{group.label}H-{indicated.group(2)}"
                indicated = None
                continue
            if not partners:
                raise UnsupportedStructure("a ketone on a doubly bonded atom without a partner for added hydrogen")
            partner = min(partners, key=sort_key)
            remaining.remove(group.label)
            remaining.remove(partner)
            suffix_locants[index] = f"{suffix_locants[index]}({partner}H)"
    if len(remaining) % 2:
        raise UnsupportedStructure("an odd number of hydro atoms would need indicated hydrogen")
    return sorted(remaining, key=sort_key), name


def _prefix_locant_key(sort_key):
    def key(locant):
        if locant.startswith("N"):
            rest = locant.lstrip("N").rstrip("′")
            return (-1, locant.count("′"), sort_key(rest) if rest else ())
        return (0, 0, sort_key(locant.rstrip("αβξ")))

    return key


def _join(parts):
    text = ""
    for part in parts:
        if not part:
            continue
        if text and (part[0].isdigit() or not part[0].isascii() or part[0] in "(["):
            text += "-"
        text += part
    return text


def n_roots(mol, group):
    """[(locant letter, nitrogen, substituent roots)] for the nitrogen atoms of an amide-type or amine-type group."""
    if group.cls in ("amide", "sulfonamide"):
        nitrogen = group.extra["nitrogen"]
        anchor = group.atoms[0] if group.cls == "amide" else group.root
        return [("N", nitrogen, _others(mol, nitrogen, anchor))]
    if group.cls == "hydrazide":
        first, last = group.extra["nitrogen"], group.extra["terminal"]
        return [("N", first, _others(mol, first, group.atoms[0], last)), ("N′", last, _others(mol, last, first))]
    if group.cls == "amidine":
        amino, imino = group.extra["nitrogen"], group.extra["imino"]
        return [("N", amino, _others(mol, amino, group.atoms[0])), ("N′", imino, _others(mol, imino, group.atoms[0]))]
    return [("N", group.atoms[0], list(group.extra["substituents"]))]


def _others(mol, nitrogen, *excluded):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(nitrogen).GetNeighbors() if n.GetIdx() not in excluded]


def _ester_alkyls(mol, graph, halogens, members, sort_key):
    names = {}
    for group in members:
        oxygen, alkyl = group.extra["oxygen"], group.extra["alkyl"]
        names[group.label] = name_branch(graph, alkyl, oxygen, halogens, frozenset(), mol=mol, unsaturated=True)
    distinct = {n for n, _ in names.values()}
    if len(distinct) == 1:
        name, compound = next(iter(names.values()))
        if len(members) > 1:
            return multiplying_prefix(len(members), compound=compound) + (f"({name})" if compound else name)
        return name
    by_name = {}
    for label, (name, compound) in names.items():
        by_name.setdefault(name, (compound, []))[1].append(label)
    parts = []
    for name, (compound, labels) in sorted(by_name.items()):
        word = name if len(labels) == 1 else multiplying_prefix(len(labels), compound=compound) + (f"({name})" if compound else name)
        parts.append(f"{','.join(sorted(labels, key=sort_key))}-{word}")
    return " ".join(parts)


def _acid_anion(mol, mapped, group, context):
    name = _acid_name(mol, mapped, group.root, group.extra["acyl"], context)
    return name[: -len("ic acid")] + "ate"


def _acid_name(mol, mapped, oxygen, acyl, context):
    """Name of the carboxylic acid R-CO-OH whose R-CO is `acyl` and whose hydroxy oxygen is `oxygen`."""
    from .core import smiles_to_iupac

    keep, stack = {acyl}, [acyl]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetIdx() != oxygen and n.GetIdx() not in keep:
                if n.GetIdx() in mapped:
                    raise UnsupportedStructure("a ring-closing O-acyl group is not supported")
                keep.add(n.GetIdx())
                stack.append(n.GetIdx())
    for atom_index in keep:
        context["used"].add(("atom", atom_index))
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in keep and bond.GetEndAtomIdx() in keep:
            context["used"].add(("bond", (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
    editable = Chem.RWMol(mol)
    ester_oxygen = editable.GetAtomWithIdx(oxygen)
    ester_oxygen.SetNumExplicitHs(1)
    ester_oxygen.SetNoImplicit(True)
    for idx in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in keep | {oxygen}), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    name = smiles_to_iupac(Chem.MolToSmiles(acid))
    if not name.endswith("ic acid"):
        raise UnsupportedStructure("the acid part of this ester is not a plain carboxylic acid")
    return name


def assemble(mol, skeleton, choice, stereo, sort_key):
    context = {"atoms": {}, "bonds": {}, "used": set()}
    for kind, index, code in stereo.outside:
        if kind == "atom":
            context["atoms"][index] = code
        else:
            bond = mol.GetBondWithIdx(index)
            context["bonds"][(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())] = code
    token = BRANCH_STEREO.set(context if stereo.outside else None)
    try:
        text = _assemble(mol, skeleton, choice, stereo, sort_key, context)
    finally:
        BRANCH_STEREO.reset(token)
    for kind, index, _ in stereo.outside:
        key = ("atom", index) if kind == "atom" else ("bond", (mol.GetBondWithIdx(index).GetBeginAtomIdx(), mol.GetBondWithIdx(index).GetEndAtomIdx()))
        if key not in context["used"]:
            raise UnsupportedStructure("a stereo element outside the parent is not cited by any substituent name")
    return text


def _assemble(mol, skeleton, choice, stereo, sort_key, context):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    mapping = choice.mapping
    mapped = set(mapping.values())
    faces = stereo.faces
    principal = None if choice.attach else principal_class(choice.groups)
    prefix_groups = {}

    def add(name, compound, locant):
        prefix_groups.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)

    members = []
    for group in choice.groups:
        if principal is not None and group.cls == principal:
            members.append(group)
            continue
        if group.kind == "o":
            raise UnsupportedStructure("a terminal functional carbon of the parent below the principal group")
        root = group.root
        name, compound = name_branch(graph, root, mapping[group.label], halogens, frozenset(), mol=mol, unsaturated=True)
        add(name, compound, _locant_text(group.label, faces, root))
    for label, root in choice.branches:
        name, compound = name_branch(graph, root, mapping[label], halogens, frozenset(), mol=mol, unsaturated=True)
        add(name, compound, _locant_text(label, faces, root))

    suffix_word, suffix_locants, tail, alkyl_front, anhydride_other = "", [], "", "", ""
    if members:
        kinds = {g.kind for g in members}
        if len(kinds) != 1:
            raise UnsupportedStructure("mixed ring and chain groups of one class are not supported")
        kind = next(iter(kinds))
        order = sorted(members, key=lambda g: sort_key(g.label))
        suffix_locants = [_locant_text(g.label, faces, g.root) for g in order]
        if principal == "ester_o":
            anions = {g.label: _acid_anion(mol, mapped, g, context) for g in order}
            distinct = set(anions.values())
            suffix_word = multiplied_word(len(order), "yl")
            if len(distinct) == 1:
                (anion,) = distinct
                tail = (multiplying_prefix(len(order), compound=False) + anion) if len(order) > 1 else anion
            else:
                tail = " ".join(f"{label}-{anion}" for label, anion in sorted(anions.items(), key=lambda kv: (kv[1], sort_key(kv[0]))))
        else:
            suffix_word = multiplied_word(len(order), SUFFIX[(principal, kind)])
            if principal == "ester":
                alkyl_front = _ester_alkyls(mol, graph, halogens, order, sort_key)
            if principal == "anhydride":
                if len(order) > 1:
                    raise UnsupportedStructure("several anhydride groups on a parent are not supported")
                g = order[0]
                other = _acid_name(mol, mapped, g.extra["oxygen"], g.extra["other"], context)
                anhydride_other = other[: -len(" acid")]
            if principal == "halide":
                words = sorted({g.extra["halide"] for g in order})
                tail = " ".join(multiplied_word(sum(g.extra["halide"] == w for g in order), w) for w in words)
            if principal in N_CLASSES:
                for g in order:
                    for letter, nitrogen, roots in n_roots(mol, g):
                        locant = letter if len(order) == 1 else f"{letter}{g.label}"
                        for root in roots:
                            name, compound = name_branch(graph, root, nitrogen, halogens, frozenset(), mol=mol, unsaturated=True)
                            add(name, compound, locant)
    if choice.attach is not None:
        suffix_word = multiplied_word(len(choice.attach), "yl")
        suffix_locants = [_locant_text(label, faces, dummy) for label, dummy in choice.attach]

    name = skeleton.name
    ending = next((e for e in ENDINGS if name.endswith(e)), None)
    use_endings = (
        ending is not None
        and not choice.lost
        and not choice.hydro
        and all(q == 1 for _, q, _ in choice.gained)
        and all(l not in choice.unsaturated_labels for pair, _, _ in choice.gained for l in pair)
        and not (ending == "anine" and any(p == 3 for _, _, p in choice.gained))
    )
    adds = []
    if not use_endings:
        dehydro = sorted((l for pair, q, p in choice.gained for _ in range(p - q) for l in pair), key=sort_key)
        hydro = sorted([l for pair in choice.lost for l in pair] + list(choice.hydro), key=sort_key)
        hydro, name = _without_added_hydrogen(hydro, skeleton, choice, members, principal, suffix_locants, sort_key)
        shared = set(hydro) & set(dehydro)
        if shared:
            left_hydro = [l for l in hydro if l not in shared]
            left_dehydro = [l for l in dehydro if l not in shared]
            if not len(left_hydro) % 2 and not len(left_dehydro) % 2:
                hydro, dehydro = left_hydro, left_dehydro
        if dehydro:
            adds.append(f"{','.join(dehydro)}-{multiplied_word(len(dehydro), 'dehydro')}")
        if hydro:
            adds.append(f"{','.join(hydro)}-{multiplied_word(len(hydro), 'hydro')}")
    ene = [_ene_text(pair, sort_key) for pair, q, p in choice.gained if use_endings and p == 2]
    yne = [_ene_text(pair, sort_key) for pair, q, p in choice.gained if use_endings and p == 3]
    if use_endings and (ene or yne):
        core = _words(name[: -len(ending)], ene, yne, ending, suffix_word, suffix_locants)
    elif suffix_word:
        core = f"{_elide(name, suffix_word)}-{','.join(suffix_locants)}-{suffix_word}"
    else:
        core = name
    core = stereo.parent + core

    if prefix_groups:
        for group in prefix_groups.values():
            group["locants"].sort(key=_prefix_locant_key(sort_key))
    prefix = format_substituent_prefixes(prefix_groups) if prefix_groups else ""
    body = stereo.front + _join([prefix] + adds + [core])
    if anhydride_other:
        first, second = sorted([anhydride_other, body], key=lambda w: re.sub(r"[^a-zα-ω]", "", w.lower()))
        return f"{first} {second} anhydride"
    text = f"{body} {tail}" if tail else body
    return f"{alkyl_front} {text}" if alkyl_front else text
