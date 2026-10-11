"""Nuclides on the nitrogen atoms of a diazene, azoxy or hydrazone parent hydride (P-82.2.1, P-82.5, P-82.6.1.1,
P-82.6.3.3): the two nitrogens of the parent are told apart only by locants, so every locant of the name is cited and
the descriptor stands in front of the parent name."""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    specified_double_bond_stereo,
    substituent_locant_set_and_citation,
)
from ._diazene import _azo_pairs, check_diazene_scope, has_diazene_shape, pair_substituents
from ._hydrazone import has_hydrazone_shape
from ._isotope_labels import descriptor, modification_key, split_isotopes
from ._multiplicative_text import enclose
from ._substituents import ISOTOPE_LABELS, format_substituent_prefixes


def _split(mol):
    split = split_isotopes(mol)
    if split is None:
        raise UnsupportedStructure("no nuclide to cite")
    return split[0], split[1]


def _shape_of_unlabelled(mol, has_shape):
    try:
        split = split_isotopes(mol)
    except UnsupportedStructure:
        return False
    return split is not None and has_shape(split[0])


def has_labelled_diazene(mol):
    return _shape_of_unlabelled(mol, has_diazene_shape)


def has_labelled_hydrazone(mol):
    return _shape_of_unlabelled(mol, has_hydrazone_shape)


def has_labelled_azoxy(mol):
    from ._dipolar import _azoxy

    return _shape_of_unlabelled(mol, lambda clean: _azoxy(clean, adjacency(clean)) is not None)


def _named_with_labels(clean, labels, parent, build):
    """`build()` run with the nuclides outside `parent` available to the substituent namers, and the nuclides of `parent`
    returned; every other nuclide must have been cited."""
    parent_labels = {a: e for a, e in labels.items() if a in parent}
    if not parent_labels:
        raise UnsupportedStructure("no nuclide on the nitrogen atoms of the parent")
    context = {"labels": {a: e for a, e in labels.items() if a not in parent}, "consumed": set(), "mol": clean}
    token = ISOTOPE_LABELS.set(context if context["labels"] else None)
    try:
        built = build()
    finally:
        ISOTOPE_LABELS.reset(token)
    if set(context["labels"]) - context["consumed"]:
        raise UnsupportedStructure("a labelled atom outside the parent is not cited by this name")
    return built, parent_labels


def _locant_text(clean, parent_labels, locant_of):
    capacity = {a: clean.GetAtomWithIdx(a).GetTotalNumHs() for a in parent_labels}
    return descriptor(parent_labels, locant_of, False, capacity=capacity)


def _orientation(substituents, parent_labels, first, second):
    """({nitrogen: locant}, grouped substituents) numbering the nitrogens 1 and 2: lowest locants to the substituents, then
    alphanumerical order, then the modified atoms (P-82.5.1, P-82.5.2)."""
    best = None
    for one, two in ((first, second), (second, first)):
        locant_of = {one: 1, two: 2}
        grouped = group_substituents({locant_of[n]: [entry] for n, entry in substituents.items()})
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        key = (locant_set, citation, modification_key(parent_labels, locant_of))
        if best is None or key < best[0]:
            best = (key, locant_of, grouped)
    return best[1], best[2]


def name_labelled_diazene(mol):
    clean, labels = _split(mol)
    pairs = _azo_pairs(clean)
    if len(pairs) != 1:
        raise UnsupportedStructure("only a single N=N pair is named with nuclides on its nitrogens")
    aromatic = check_diazene_scope(clean, pairs)
    first, second = (a.GetIdx() for a in pairs[0])
    graph = adjacency(clean)
    substituents, parent_labels = _named_with_labels(
        clean, labels, {first, second},
        lambda: pair_substituents(clean, graph, halogen_substituents(clean), aromatic, first, second),
    )
    locant_of, grouped = _orientation(substituents, parent_labels, first, second)
    return format_substituent_prefixes(grouped) + _locant_text(clean, parent_labels, locant_of) + "diazene"


def name_labelled_azoxy(mol):
    from ._dipolar import _azoxy, _group

    clean, labels = _split(mol)
    graph = adjacency(clean)
    unit = _azoxy(clean, graph)
    if (
        unit is None
        or len(Chem.GetMolFrags(clean)) != 1
        or sum(1 for a in clean.GetAtoms() if a.GetFormalCharge()) != 2
        or any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in clean.GetBonds())
        or any(a.GetAtomicNum() not in (6, 7, 8, 9, 17, 35, 53) for a in clean.GetAtoms())
    ):
        raise UnsupportedStructure("this azoxy compound is not named with nuclides on its nitrogens")
    oxide_nitrogen, other_nitrogen, oxygen = (a.GetIdx() for a in unit)
    if oxygen in labels:
        raise UnsupportedStructure("a labelled oxide oxygen is not supported yet")

    def groups():
        found = {}
        for nitrogen in (oxide_nitrogen, other_nitrogen):
            (root,) = [n for n in graph[nitrogen] if n not in (oxygen, oxide_nitrogen, other_nitrogen)]
            found[nitrogen] = _group(clean, graph, root, nitrogen)
        return found

    substituents, parent_labels = _named_with_labels(clean, labels, {oxide_nitrogen, other_nitrogen}, groups)
    if substituents[oxide_nitrogen] == substituents[other_nitrogen]:
        locant_of = min(
            ({oxide_nitrogen: 1, other_nitrogen: 2}, {other_nitrogen: 1, oxide_nitrogen: 2}),
            key=lambda option: (option[oxide_nitrogen], modification_key(parent_labels, option)),
        )
        grouped = group_substituents({locant_of[n]: [entry] for n, entry in substituents.items()})
    else:
        locant_of, grouped = _orientation(substituents, parent_labels, oxide_nitrogen, other_nitrogen)
    prefixes = format_substituent_prefixes(grouped)
    return f"{prefixes}{_locant_text(clean, parent_labels, locant_of)}diazene {locant_of[oxide_nitrogen]}-oxide"


def name_labelled_hydrazone(mol):
    from ._hydrazone import _name_hydrazone_carbon, _validate_and_find_hydrazone

    clean, labels = _split(mol)
    carbon, imine_nitrogen = _validate_and_find_hydrazone(clean)
    (amino,) = [n.GetIdx() for n in clean.GetAtomWithIdx(imine_nitrogen).GetNeighbors() if n.GetIdx() != carbon]
    name, parent_labels = _named_with_labels(
        clean, labels, {imine_nitrogen, amino}, lambda: _name_hydrazone_carbon(clean, carbon, imine_nitrogen)
    )
    group = enclose(name) if any(ch.isdigit() or ch in "(-" for ch in name) else name
    text = f"1-{group}{_locant_text(clean, parent_labels, {imine_nitrogen: 1, amino: 2})}hydrazine"
    stereo = specified_double_bond_stereo(clean)
    return f"({stereo[0][1]})-{text}" if stereo is not None else text
