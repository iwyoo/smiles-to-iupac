"""Secondary and tertiary amides (R-CO)2NH, (R-SO2)2NH, (R-CO)3N (P-66.1.2): N-acyl derivatives of the senior primary
amide, e.g. N-acetylbenzamide, N,N-di(cyclohexanecarbonyl)cyclohexanecarboxamide. The other acyl groups are stood in
for by methyl groups whose names are replaced by the acyl prefixes, so the amide is named with its usual N-prefixes."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._diester_ring_diyl import _system_of
from ._dipolar import _smiles_with_order
from ._ring_system_seniority import ring_seniority_key
from ._substituents import FORCED_BRANCH_NAMES, name_branch


def _acyl_roots(mol, nitrogen):
    """Neighbours of `nitrogen` that are the carbon of a carboxylic acyl group or the sulfur of a sulfonyl group."""
    roots = []
    for n in nitrogen.GetNeighbors():
        if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        others = [m for m in n.GetNeighbors() if m.GetIdx() != nitrogen.GetIdx()]
        doubly = [m for m in others if m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 2.0]
        rest = [m for m in others if m not in doubly]
        if n.GetAtomicNum() == 6 and len(doubly) == 1 and not n.IsInRing():
            if (not rest and n.GetTotalNumHs() == 1) or (len(rest) == 1 and rest[0].GetAtomicNum() == 6):
                roots.append(n.GetIdx())
        elif n.GetAtomicNum() == 16 and len(doubly) == 2 and len(rest) == 1 and rest[0].GetAtomicNum() == 6 and not n.IsInRing():
            roots.append(n.GetIdx())
    return roots


def _side(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _rank(mol, graph, root, nitrogen):
    """Smaller is senior: carboxamides before sulfonamides, ring parents before chains, heterocycles before
    carbocycles, then the larger ring or the longer chain."""
    atom = mol.GetAtomWithIdx(root)
    sulfur = atom.GetAtomicNum() == 16
    carbons = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    if not carbons:
        return (int(sulfur), 2, (), 0)
    start = carbons[0]
    side = _side(graph, start, root)
    if mol.GetAtomWithIdx(start).IsInRing():
        return (int(sulfur), 0, ring_seniority_key(mol, _system_of(mol, start)[1]), 0)
    length = len([a for a in side if mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(a).IsInRing()])
    return (int(sulfur), 1, (), -length)


def diacylamine_name(mol):
    """The name of a molecule with an acyclic nitrogen bearing two or three acyl groups, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return None
    nitrogens = [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 7 and not a.IsInRing() and not a.GetIsAromatic() and len(_acyl_roots(mol, a)) >= 2
    ]
    if len(nitrogens) != 1:
        return None
    nitrogen = nitrogens[0]
    graph = adjacency(mol)
    roots = _acyl_roots(mol, nitrogen)
    parent = min(roots, key=lambda r: (_rank(mol, graph, r, nitrogen.GetIdx()), r))
    donors = [r for r in roots if r != parent]
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    try:
        names = {r: name_branch(graph, r, nitrogen.GetIdx(), halogens, aromatic, mol=mol) for r in donors}
    except UnsupportedStructure:
        return None
    editable = Chem.RWMol(mol)
    placeholders = {}
    removed = set()
    for r in donors:
        side = _side(graph, r, nitrogen.GetIdx())
        removed |= side - {r}
        editable.RemoveBond(nitrogen.GetIdx(), r)
    keep = sorted(removed, reverse=True)
    for donor in donors:
        atom = editable.GetAtomWithIdx(donor)
        atom.SetAtomicNum(6)
        for bond in list(atom.GetBonds()):
            if bond.GetOtherAtomIdx(donor) in removed:
                editable.RemoveBond(donor, bond.GetOtherAtomIdx(donor))
        editable.AddBond(nitrogen.GetIdx(), donor, Chem.BondType.SINGLE)
        placeholders[donor] = None
    for idx in keep:
        editable.RemoveAtom(idx)
    mapping = {}
    shift = 0
    old_to_new = {}
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
    forced = {position[old_to_new[r]]: names[r] for r in donors}
    from .core import smiles_to_iupac

    token = FORCED_BRANCH_NAMES.set((stand_in.GetNumAtoms(), forced))
    try:
        name = smiles_to_iupac(smiles)
    except UnsupportedStructure:
        return None
    finally:
        FORCED_BRANCH_NAMES.reset(token)
    if not name.endswith("amide") or any(f"{donor_name}" not in name for donor_name, _ in names.values()):
        return None
    return name
