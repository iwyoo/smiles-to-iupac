"""Amine and imine oxides and their chalcogen analogues (P-62.5): functional class nomenclature for one oxide,
'N,N-dimethylmethanamine N-oxide'; an oxide on a further nitrogen is cited as an '(oxo-λ5-azanyl)' prefix."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents, with_unlocated_descriptor
from ._phosphanyl_group import PREFIX_PROP
from ._prefix_groups import enclose
from ._substituents import name_branch


_CHALCOGEN_CLASS = {8: "oxide", 16: "sulfide", 34: "selenide", 52: "telluride"}
_NUCLIDE_GROUP = re.compile(r"\((?:[^()]*-)?(?:\d+[A-Z][a-z]?\d*,?)+\)")


def _unlabelled_text(name):
    return re.sub(r"[()\[\]{}]", "", _NUCLIDE_GROUP.sub("", name))


def _oxide_ligand(mol, nitrogen):
    """The terminal chalcogenide atom of an amine or imine oxide nitrogen, else None."""
    if nitrogen.GetAtomicNum() != 7 or nitrogen.GetFormalCharge() != 1 or nitrogen.GetIsotope() != 0:
        return None
    if nitrogen.GetDegree() + nitrogen.GetTotalNumHs() not in (3, 4):
        return None
    ligands = [
        n
        for n in nitrogen.GetNeighbors()
        if n.GetAtomicNum() in _CHALCOGEN_CLASS and n.GetFormalCharge() == -1 and n.GetDegree() == 1
    ]
    if len(ligands) != 1 or ligands[0].GetIsotope() != 0:
        return None
    (ligand,) = ligands
    if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), ligand.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    others = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != ligand.GetIdx()]
    if len(others) not in (1, 2, 3) or any(n.GetAtomicNum() != 6 for n in others):
        return None
    if any(mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in others):
        return None
    return ligand


def _oxide_nitrogens(mol):
    return [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7 and a.GetFormalCharge() == 1]


def has_amine_oxide_shape(mol) -> bool:
    nitrogens = _oxide_nitrogens(mol)
    return bool(nitrogens) and all(_oxide_ligand(mol, n) is not None for n in nitrogens)


def name_amine_oxide(mol) -> str:
    from ._polyfunctional import name_polyfunctional

    nitrogens = _oxide_nitrogens(mol)
    if len(nitrogens) > 1:
        return _name_among_oxides(mol, nitrogens)
    nitrogen = nitrogens[0]
    oxide_oxygen = _oxide_ligand(mol, nitrogen)
    class_word = _CHALCOGEN_CLASS[oxide_oxygen.GetAtomicNum()]

    reduced_mol = _reduce(mol, nitrogen, oxide_oxygen)
    if not _is_amine_nitrogen(reduced_mol.GetAtomWithIdx(_reduced_position(nitrogen.GetIdx(), oxide_oxygen.GetIdx()))):
        raise UnsupportedStructure("the oxidized nitrogen carries an acyl or nitrile carbon, whose oxide is not named here")

    amine_nitrogens = [a for a in reduced_mol.GetAtoms() if _is_amine_nitrogen(a)]
    other_nitrogens = sum(a.GetAtomicNum() == 7 for a in reduced_mol.GetAtoms()) - len(amine_nitrogens)
    if len(amine_nitrogens) == 1:
        labelled = _labelled_oxide_name(mol, nitrogen, oxide_oxygen, class_word)
        if labelled is not None:
            return labelled
    if len(amine_nitrogens) == 1 and other_nitrogens:
        if nitrogen.IsInRing():
            raise UnsupportedStructure("the ring of an oxidized ring nitrogen is the parent of its N-oxide (P-62.5)")
        return f"{_name_with_nitrile_prefix(reduced_mol)} N-{class_word}"
    if len(amine_nitrogens) > 1:
        if nitrogen.IsInRing():
            raise UnsupportedStructure("the ring of an oxidized ring nitrogen is the parent of its N-oxide (P-62.5)")
        oxide_idx = oxide_oxygen.GetIdx()
        named = _in_reduced_molecule(
            mol, reduced_mol, oxide_idx, lambda reduced: _name_with_oxidized_parent(reduced, nitrogen.GetIdx(), oxide_idx)
        )
        return f"{named} N-{class_word}"
    from ._common import match_cip_labels, specified_stereo_elements

    nitrogen_code = next(
        (code for kind, idx, code in specified_stereo_elements(mol) or () if kind == "atom" and idx == nitrogen.GetIdx()),
        None,
    )
    reduced_mol.GetAtomWithIdx(_reduced_position(nitrogen.GetIdx(), oxide_oxygen.GetIdx())).SetChiralTag(
        Chem.ChiralType.CHI_UNSPECIFIED
    )
    match_cip_labels(
        mol,
        reduced_mol,
        {i: i + (i >= oxide_oxygen.GetIdx()) for i in range(reduced_mol.GetNumAtoms())},
    )
    if nitrogen.IsInRing():
        raise UnsupportedStructure("the ring of an oxidized ring nitrogen is the parent of its N-oxide (P-62.5)")
    try:
        base_name = name_polyfunctional(reduced_mol)
    except UnsupportedStructure:
        base_name = _name_with_nitrile_prefix(reduced_mol)
    return f"{with_unlocated_descriptor(base_name, nitrogen_code)} N-{class_word}"


def _in_reduced_molecule(mol, reduced_mol, oxide_idx, name_it):
    """`name_it(reduced_mol)` with the nuclide labels of `mol` moved to the atoms of `reduced_mol`, which lacks the oxide
    atom, and the atoms it cites marked in the labels of `mol`."""
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if not context or not context["labels"] or context["mol"] is not mol:
        return name_it(reduced_mol)
    if oxide_idx in context["labels"]:
        raise UnsupportedStructure("a nuclide on the oxide group of an amine oxide is not supported yet")

    def shifted(idx):
        return idx - (idx > oxide_idx)

    inner = {"labels": {shifted(i): entry for i, entry in context["labels"].items()}, "consumed": set(), "mol": reduced_mol}
    token = ISOTOPE_LABELS.set(inner)
    try:
        name = name_it(reduced_mol)
    finally:
        ISOTOPE_LABELS.reset(token)
    original = {shifted(i): i for i in context["labels"]}
    context["consumed"].update(original[i] for i in inner["consumed"] if i in original)
    return name


def _labelled_oxide_name(mol, nitrogen, oxide_oxygen, class_word):
    """The name of the amine with the nuclides of `mol` followed by the class word, or None without nuclides; the nuclides
    go in front of the amine name as for the neutral amine (P-82.2.1, P-62.5)."""
    from ._isotope_labels import restore_labels
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if not context or not context["labels"] or context["mol"] is not mol:
        return None
    labels = context["labels"]
    oxide = labels.get(oxide_oxygen.GetIdx())
    editable = Chem.RWMol(mol)
    restore_labels(editable, mol, labels, skip={oxide_oxygen.GetIdx()})
    editable.GetAtomWithIdx(nitrogen.GetIdx()).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    editable.GetAtomWithIdx(nitrogen.GetIdx()).SetFormalCharge(0)
    editable.RemoveAtom(oxide_oxygen.GetIdx())
    reduced = editable.GetMol()
    Chem.SanitizeMol(reduced)
    token = ISOTOPE_LABELS.set(None)
    try:
        name = _name_with_nitrile_prefix(reduced)
    finally:
        ISOTOPE_LABELS.reset(token)
    context["consumed"].update(labels)
    return f"{name} N-{'(' + oxide['skeleton'] + ')' if oxide else ''}{class_word}"


def _is_amine_nitrogen(atom):
    """A nitrogen with single bonds to carbon or hydrogen only, none of the carbons an acyl or nitrile carbon."""
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge():
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondTypeAsDouble() != 1.0 or other.GetAtomicNum() != 6:
            return False
        if any(b.GetBondTypeAsDouble() > 1.0 and b.GetOtherAtom(other).GetAtomicNum() != 6 for b in other.GetBonds()):
            return False
    return True


def _name_with_nitrile_prefix(reduced_mol):
    """P-62.5, P-67.1.6: the amine oxide is the parent, so a nitrile is cited as the prefix 'cyano' ('cyano-N,N-
    dimethylmethanamine N-oxide') and a carboxy group beside it as 'carboxy'."""
    from ._polyfunctional import FORCED_PRINCIPAL, name_polyfunctional

    token = FORCED_PRINCIPAL.set("amine")
    try:
        name = name_polyfunctional(reduced_mol)
    finally:
        FORCED_PRINCIPAL.reset(token)
    # P-14.3.4.4: the oxidized nitrogen carries no hydrogen, so the only position left for the prefix is carbon 1
    cited_prefix = any(a.HasProp(PREFIX_PROP) for a in reduced_mol.GetAtoms())
    return re.sub(r"^1-(?=[a-z(\[{])", "", name) if name.endswith("methanamine") and not cited_prefix else name


def _reduce(mol, nitrogen, oxide_oxygen):
    reduced = Chem.RWMol(mol)
    reduced.GetAtomWithIdx(nitrogen.GetIdx()).SetFormalCharge(0)
    reduced.RemoveAtom(oxide_oxygen.GetIdx())
    reduced_mol = reduced.GetMol()
    Chem.SanitizeMol(reduced_mol)
    return reduced_mol


def _reduced_position(nitrogen_idx, oxide_idx):
    return nitrogen_idx - (1 if oxide_idx < nitrogen_idx else 0)


def _name_among_oxides(mol, nitrogens):
    """P-62.5: one oxide gives the class term, each further one is an '(oxo-λ5-azanyl)' prefix; the parent amine is the
    one whose carbon skeleton is senior (ring before chain, then the larger acyclic carbon set, then the senior chalcogen),
    of equal ones the more modified (P-44.4.1.11.1)."""
    from ._polyfunctional import name_polyfunctional
    from ._substituents import ISOTOPE_LABELS
    from .core import _name_mol

    graph = adjacency(mol)
    outer = ISOTOPE_LABELS.get()
    labels = outer["labels"] if outer and outer["mol"] is mol else {}
    start = set(outer["consumed"]) if labels else set()
    ranked = []
    for parent in nitrogens:
        if labels:
            outer["consumed"].clear()
            outer["consumed"].update(start)
        found = _cite_other_oxides(mol, graph, parent, nitrogens)
        if found is None:
            continue
        contracted, removed = found
        kept = {i: entry for i, entry in labels.items() if i not in removed}
        if labels and any(i not in outer["consumed"] for i in labels if i in removed):
            continue
        shift = {i: i - sum(1 for r in removed if r < i) for i in kept}
        inner = {"labels": {shift[i]: entry for i, entry in kept.items()}, "consumed": set(), "mol": contracted}
        token = ISOTOPE_LABELS.set(inner if kept else None)
        try:
            name = _name_mol(contracted)
        except UnsupportedStructure:
            name = None
            if parent.IsInRing():
                try:
                    name = name_polyfunctional(contracted)
                except UnsupportedStructure:
                    pass
        finally:
            ISOTOPE_LABELS.reset(token)
        if name is None or set(inner["labels"]) - inner["consumed"]:
            continue
        ranked.append(((*_skeleton_rank(mol, graph, parent), len(kept)), name, set(outer["consumed"]) | set(kept) if labels else set()))
    if not ranked:
        raise UnsupportedStructure("no amine parent carries one of the oxidized nitrogens of this polyamine oxide")
    best = max(rank for rank, _, _ in ranked)
    name, consumed = min((name, consumed) for rank, name, consumed in ranked if rank == best)
    if labels:
        outer["consumed"].clear()
        outer["consumed"].update(consumed)
    return name


def _cite_other_oxides(mol, graph, parent, nitrogens):
    """(`mol` with every oxidized nitrogen but `parent` collapsed to a placeholder carrying its prefix, the atoms removed
    for them), or None."""
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    groups, removed = [], set()
    for other in nitrogens:
        idx = other.GetIdx()
        if idx == parent.GetIdx() or idx in removed:
            continue
        toward = [n for n in graph[idx] if n in _reachable(graph, parent.GetIdx(), {idx})]
        if len(toward) != 1 or other.IsInRing():
            return None
        try:
            name, compound = name_branch(graph, idx, toward[0], halogens, aromatic, mol=mol)
        except UnsupportedStructure:
            return None
        atoms = _reachable(graph, idx, {toward[0]})
        if parent.GetIdx() in atoms:
            return None
        groups.append((toward[0], enclose(name) if compound else name))
        removed |= atoms
    rw = Chem.RWMol(mol)
    for anchor, name in groups:
        placeholder = rw.AddAtom(Chem.Atom(53))
        rw.AddBond(anchor, placeholder, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(placeholder).SetProp(PREFIX_PROP, name)
    for idx in sorted(removed, reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out, removed


def _reachable(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n not in blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _skeleton_rank(mol, graph, parent):
    def acyclic_carbon(idx):
        atom = mol.GetAtomWithIdx(idx)
        return atom.GetAtomicNum() == 6 and not atom.IsInRing()

    carbons = [n for n in graph[parent.GetIdx()] if mol.GetAtomWithIdx(n).GetAtomicNum() == 6]
    blocked = {i for i in range(mol.GetNumAtoms()) if not acyclic_carbon(i)}
    longest = max((len(_reachable(graph, c, blocked)) for c in carbons if acyclic_carbon(c)), default=0)
    ring = any(mol.GetAtomWithIdx(c).IsInRing() for c in carbons)
    return ring, longest, -_oxide_ligand(mol, parent).GetAtomicNum()


def _name_with_oxidized_parent(reduced_mol, nitrogen_idx, oxide_idx):
    """P-62.5: the oxidized nitrogen is the amine suffix nitrogen of the parent, so every other amino group is cited
    as a prefix ('5-(dimethylamino)-N,N-dimethylpentan-1-amine N-oxide')."""
    from .core import _name_mol, nuclides_cited
    from ._hetero_chain import KEPT_ATOMS_PROP, contract_hetero_groups_candidates
    from ._isotope_labels import restore_labels
    from ._polyfunctional import name_polyfunctional
    from ._substituents import ISOTOPE_LABELS

    position = _reduced_position(nitrogen_idx, oxide_idx)
    outer = ISOTOPE_LABELS.get()
    labels = outer["labels"] if outer and outer["mol"] is reduced_mol else {}
    start = set(outer["consumed"]) if labels else set()
    names = []
    try:
        candidates = contract_hetero_groups_candidates(reduced_mol, {position})
    except UnsupportedStructure:
        if not labels:
            raise
        candidates = []
    for contracted in candidates:
        if not labels:
            try:
                names.append((_name_mol(contracted), set()))
            except UnsupportedStructure:
                pass
            continue
        kept = [int(i) for i in contracted.GetProp(KEPT_ATOMS_PROP).split(",")]
        new_of = {old: new for new, old in enumerate(kept)}
        carried = {i: entry for i, entry in labels.items() if i in new_of}
        if any(i not in outer["consumed"] for i in labels if i not in new_of):
            continue
        editable = Chem.RWMol(contracted)
        restore_labels(editable, contracted, {new_of[i]: entry for i, entry in carried.items()})
        labelled = editable.GetMol()
        Chem.SanitizeMol(labelled)
        token = ISOTOPE_LABELS.set(None)
        try:
            name = _name_mol(labelled)
        except UnsupportedStructure:
            continue
        finally:
            ISOTOPE_LABELS.reset(token)
        if not nuclides_cited(labelled, name):
            continue
        names.append((name, set(outer["consumed"]) | set(carried)))
    if labels:
        outer["consumed"].clear()
        outer["consumed"].update(start)
    if not names:
        reduced_mol.GetAtomWithIdx(position).SetBoolProp("_oxidized_amine", True)
        target = reduced_mol
        if labels:
            editable = Chem.RWMol(reduced_mol)
            restore_labels(editable, reduced_mol, labels)
            target = editable.GetMol()
            Chem.SanitizeMol(target)
        token = ISOTOPE_LABELS.set(None)
        try:
            name = name_polyfunctional(target)
            if labels:
                if not nuclides_cited(target, name):
                    raise UnsupportedStructure("the nuclides of this polyamine N-oxide are not cited by the name of its parent")
                if _unlabelled_text(name) != _unlabelled_text(name_polyfunctional(reduced_mol)):
                    raise UnsupportedStructure("the nuclides choose another parent than the oxidized nitrogen")
                outer["consumed"].update(labels)
            return name
        except UnsupportedStructure:
            raise UnsupportedStructure("no amine parent carries the oxidized nitrogen of this polyamine N-oxide")
        finally:
            ISOTOPE_LABELS.reset(token)
    name, consumed = min(names)
    if labels:
        outer["consumed"].update(consumed)
    return name
