"""Amine and imine oxides and their chalcogen analogues (P-62.5): functional class nomenclature for one oxide,
'N,N-dimethylmethanamine N-oxide'; an oxide on a further nitrogen is cited as an '(oxo-λ5-azanyl)' prefix."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._phosphanyl_group import PREFIX_PROP
from ._prefix_groups import enclose
from ._substituents import name_branch


_CHALCOGEN_CLASS = {8: "oxide", 16: "sulfide", 34: "selenide", 52: "telluride"}


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
    from ._amine import name_amine

    nitrogens = _oxide_nitrogens(mol)
    if len(nitrogens) > 1:
        return _name_among_oxides(mol, nitrogens)
    nitrogen = nitrogens[0]
    oxide_oxygen = _oxide_ligand(mol, nitrogen)
    class_word = _CHALCOGEN_CLASS[oxide_oxygen.GetAtomicNum()]

    reduced_mol = _reduce(mol, nitrogen, oxide_oxygen)

    amine_nitrogens = [a for a in reduced_mol.GetAtoms() if _is_amine_nitrogen(a)]
    other_nitrogens = sum(a.GetAtomicNum() == 7 for a in reduced_mol.GetAtoms()) - len(amine_nitrogens)
    if len(amine_nitrogens) == 1 and other_nitrogens:
        return f"{_name_with_nitrile_prefix(reduced_mol)} N-{class_word}"
    if len(amine_nitrogens) > 1:
        return f"{_name_with_oxidized_parent(reduced_mol, nitrogen.GetIdx(), oxide_oxygen.GetIdx())} N-{class_word}"
    from ._common import heteroatom_stereo_prefix, specified_stereocenters

    centres = specified_stereocenters(mol) or []
    stereo_prefix = heteroatom_stereo_prefix(mol, nitrogen.GetIdx()) or "" if any(i == nitrogen.GetIdx() for i, _ in centres) else ""
    try:
        base_name = name_amine(reduced_mol)
    except UnsupportedStructure:
        if nitrogen.IsInRing():
            raise
        base_name = _name_with_nitrile_prefix(reduced_mol)
    return f"{stereo_prefix}{base_name} N-{class_word}"


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
    return re.sub(r"^1-(?=[a-z(\[{])", "", name) if name.endswith("methanamine") else name


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
    one whose carbon skeleton is senior (ring before chain, then the larger acyclic carbon set, then the senior chalcogen)."""
    from .core import _name_mol

    graph = adjacency(mol)
    ranked = []
    for parent in nitrogens:
        contracted = _cite_other_oxides(mol, graph, parent, nitrogens)
        if contracted is None:
            continue
        try:
            name = _name_mol(contracted)
        except UnsupportedStructure:
            continue
        ranked.append((_skeleton_rank(mol, graph, parent), name))
    if not ranked:
        raise UnsupportedStructure("no amine parent carries one of the oxidized nitrogens of this polyamine oxide")
    best = max(rank for rank, _ in ranked)
    return min(name for rank, name in ranked if rank == best)


def _cite_other_oxides(mol, graph, parent, nitrogens):
    """`mol` with every oxidized nitrogen but `parent` collapsed to a placeholder carrying its prefix, or None."""
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
    return out


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
    from .core import _name_mol
    from ._hetero_chain import contract_hetero_groups_candidates
    from ._polyfunctional import name_polyfunctional

    position = _reduced_position(nitrogen_idx, oxide_idx)
    names = []
    for contracted in contract_hetero_groups_candidates(reduced_mol, {position}):
        try:
            names.append(_name_mol(contracted))
        except UnsupportedStructure:
            continue
    if not names:
        reduced_mol.GetAtomWithIdx(position).SetBoolProp("_oxidized_amine", True)
        try:
            return name_polyfunctional(reduced_mol)
        except UnsupportedStructure:
            raise UnsupportedStructure("no amine parent carries the oxidized nitrogen of this polyamine N-oxide")
    return min(names)
