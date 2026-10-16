"""Diacyl and triacyl derivatives of hydrazine (P-66.3.3.3): N'-benzoylbenzohydrazide, N'-acetyl-N'-ethyl-N-methyl-
propanehydrazide, N'1,N'4-diacetylbutanedihydrazide. The senior acyl group is the hydrazide; the other acyl groups are
stood in for by methyl groups whose names are replaced by the acyl prefixes, as for the diacylamines."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._diacylamine import _acyl_roots, _rank, _side, donor_names, forced_for_unlabelled, modified_atoms_beyond, unlabelled_with_labels
from ._dipolar import _smiles_with_order
from ._isotope_labels import restore_labels
from ._substituents import FORCED_BRANCH_NAMES


def _hydrazine_pair(mol):
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if (
            a.GetAtomicNum() == 7
            and b.GetAtomicNum() == 7
            and bond.GetBondTypeAsDouble() == 1.0
            and not bond.IsInRing()
            and not a.GetIsAromatic()
            and not b.GetIsAromatic()
        ):
            yield a, b


def _hydrazide_groups(mol, root):
    """The carbonyl carbons of the carbon skeleton `root` belongs to: the principal groups its hydrazide would carry."""
    seen, stack = {root}, [root]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIdx() not in seen:
                seen.add(n.GetIdx())
                stack.append(n.GetIdx())
    return sum(
        1
        for c in seen
        if any(o.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(c, o.GetIdx()).GetBondTypeAsDouble() == 2.0 for o in mol.GetAtomWithIdx(c).GetNeighbors())
        and any(n.GetAtomicNum() == 7 for n in mol.GetAtomWithIdx(c).GetNeighbors())
    )


def _parent_rank(mol, graph, root, nitrogen):
    """P-44.1.1: the acid with the most principal groups gives the parent hydrazide, then the seniority of its skeleton."""
    rank = _rank(mol, graph, root, nitrogen)
    return (rank[0], -_hydrazide_groups(mol, root)) + rank[1:]


def diacylhydrazine_name(mol):
    """The name of a molecule with an acyclic N-N bond whose nitrogens carry two or more acyl groups, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    mol, labels = unlabelled_with_labels(mol)
    if mol is None:
        return None
    pairs = []
    for a, b in _hydrazine_pair(mol):
        roots = {a.GetIdx(): _acyl_roots(mol, a), b.GetIdx(): _acyl_roots(mol, b)}
        if sum(len(r) for r in roots.values()) >= 2:
            pairs.append((a, b, roots))
    nitrogens = [n for _, _, roots in pairs for n in roots]
    if not pairs or len(set(nitrogens)) != len(nitrogens):
        return None
    graph = adjacency(mol)
    parents, donors = set(), []
    for _, _, roots in pairs:
        candidates = [(r, n) for n, rs in roots.items() for r in rs]
        parent = min(
            candidates,
            key=lambda item: (_parent_rank(mol, graph, item[0], item[1]), modified_atoms_beyond(graph, labels, *item), item[0]),
        )[0]
        parents.add(parent)
        donors += [(r, n) for r, n in candidates if r != parent]
    taken = [_side(graph, r, n) for r, n in donors]
    if sum(len(s) for s in taken) != len(set().union(*taken)) or any(
        s & (parents | set(nitrogens)) for s in taken
    ):
        return None
    halogens = halogen_substituents(mol)
    aromatic = frozenset(x.GetIdx() for x in mol.GetAtoms() if x.GetIsAromatic())
    try:
        names = donor_names(mol, graph, donors, halogens, aromatic, labels)
    except UnsupportedStructure:
        return None
    editable = Chem.RWMol(mol)
    removed = set()
    for r, n in donors:
        removed |= _side(graph, r, n) - {r}
        editable.RemoveBond(n, r)
    for r, n in donors:
        atom = editable.GetAtomWithIdx(r)
        atom.SetAtomicNum(6)
        for bond in list(atom.GetBonds()):
            if bond.GetOtherAtomIdx(r) in removed:
                editable.RemoveBond(r, bond.GetOtherAtomIdx(r))
        editable.AddBond(n, r, Chem.BondType.SINGLE)
    restore_labels(editable, mol, labels, skip=removed | {r for r, _ in donors})
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    old_to_new, shift = {}, 0
    for i in range(mol.GetNumAtoms()):
        if i in removed:
            shift += 1
        else:
            old_to_new[i] = i - shift
    stand_in = editable.GetMol()
    try:
        Chem.SanitizeMol(stand_in)
    except Chem.rdchem.MolSanitizeException:
        return None
    smiles, position = _smiles_with_order(stand_in)
    forced = forced_for_unlabelled(stand_in, position, {old_to_new[r]: names[r] for r, _ in donors})
    from .core import smiles_to_iupac

    token = FORCED_BRANCH_NAMES.set(forced)
    try:
        name = smiles_to_iupac(smiles)
    except UnsupportedStructure:
        return None
    finally:
        FORCED_BRANCH_NAMES.reset(token)
    if not name.endswith("hydrazide") or any(donor_name not in name for donor_name, _ in names.values()):
        return None
    return name
