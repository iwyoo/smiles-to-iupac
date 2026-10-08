"""Diacyl and triacyl derivatives of hydrazine (P-66.3.3.3): N'-benzoylbenzohydrazide, N'-acetyl-N'-ethyl-N-methyl-
propanehydrazide. The senior acyl group is the hydrazide; the other acyl groups are stood in for by methyl groups whose
names are replaced by the acyl prefixes, as for the diacylamines."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._diacylamine import _acyl_roots, _rank, _side
from ._dipolar import _smiles_with_order
from ._substituents import FORCED_BRANCH_NAMES, name_branch


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


def diacylhydrazine_name(mol):
    """The name of a molecule with an acyclic N-N bond whose nitrogens carry two or more acyl groups, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return None
    pairs = []
    for a, b in _hydrazine_pair(mol):
        roots = {a.GetIdx(): _acyl_roots(mol, a), b.GetIdx(): _acyl_roots(mol, b)}
        if sum(len(r) for r in roots.values()) >= 2:
            pairs.append((a, b, roots))
    if len(pairs) != 1:
        return None
    a, b, roots = pairs[0]
    graph = adjacency(mol)
    candidates = [(r, n) for n, rs in roots.items() for r in rs]
    parent, parent_n = min(candidates, key=lambda item: (_rank(mol, graph, item[0], item[1]), item[0]))
    donors = [(r, n) for r, n in candidates if r != parent]
    halogens = halogen_substituents(mol)
    aromatic = frozenset(x.GetIdx() for x in mol.GetAtoms() if x.GetIsAromatic())
    try:
        names = {r: name_branch(graph, r, n, halogens, aromatic, mol=mol) for r, n in donors}
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
    forced = {position[old_to_new[r]]: names[r] for r, _ in donors}
    from .core import smiles_to_iupac

    token = FORCED_BRANCH_NAMES.set((stand_in.GetNumAtoms(), forced))
    try:
        name = smiles_to_iupac(smiles)
    except UnsupportedStructure:
        return None
    finally:
        FORCED_BRANCH_NAMES.reset(token)
    if not name.endswith("hydrazide") or any(donor_name not in name for donor_name, _ in names.values()):
        return None
    return name
