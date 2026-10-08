"""Acylated ends of a homogeneous nitrogen chain of three or more atoms (P-58.3.2): the chain is broken so the acyl
nitrogen is expressed as an amide, 'N-(triazan-1-yl)benzamide'; two identical acylated ends joined by the rest of the
chain are named multiplicatively, 'N,N'-(hydrazine-1,2-diyl)dibenzamide'. The rest of the chain is stood in for by a
carbon whose name is replaced, as for the diacylamines."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._diacylamine import _acyl_roots, _side
from ._dipolar import _smiles_with_order
from ._numerals import multiplying_prefix
from ._substituents import FORCED_BRANCH_NAMES, name_branch

_LINKER = {1: ("methylene", "azanediyl"), 2: ("ethane-1,2-diyl", "hydrazine-1,2-diyl")}


def _plain_nitrogen(atom):
    return atom.GetAtomicNum() == 7 and not atom.IsInRing() and not atom.GetIsAromatic() and not atom.GetFormalCharge()


def _chain_partner(mol, nitrogen, acyl):
    partners = [
        n
        for n in nitrogen.GetNeighbors()
        if n.GetIdx() != acyl and _plain_nitrogen(n) and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
    ]
    return partners[0] if len(partners) == 1 else None


def _extends(mol, partner, behind):
    return any(n.GetAtomicNum() == 7 and n.GetIdx() != behind for n in partner.GetNeighbors())


def _acylated_ends(mol):
    """(nitrogen, acyl root, chain partner) for every acylated nitrogen whose N-N partner continues the chain."""
    ends = []
    for atom in mol.GetAtoms():
        if not _plain_nitrogen(atom):
            continue
        roots = _acyl_roots(mol, atom)
        if len(roots) != 1 or not any(n.GetAtomicNum() != 1 for n in mol.GetAtomWithIdx(roots[0]).GetNeighbors() if n.GetIdx() != atom.GetIdx() and n.GetAtomicNum() != 8):
            continue
        partner = _chain_partner(mol, atom, roots[0])
        if partner is not None and _extends(mol, partner, atom.GetIdx()):
            ends.append((atom.GetIdx(), roots[0], partner.GetIdx()))
    return ends


def _stand_in(mol, replacements):
    """`mol` with each (anchor, root) side replaced by one carbon on anchor, and the forced-name positions of those carbons."""
    graph = adjacency(mol)
    editable = Chem.RWMol(mol)
    removed = set()
    for anchor, root in replacements:
        removed |= _side(graph, root, anchor)
    carbons = {}
    for anchor, root in replacements:
        carbons[root] = editable.AddAtom(Chem.Atom(6))
        editable.AddBond(anchor, carbons[root], Chem.BondType.SINGLE)
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    shift = lambda i: i - sum(1 for r in removed if r < i)
    stand_in = editable.GetMol()
    Chem.SanitizeMol(stand_in)
    return stand_in, {root: shift(index) for root, index in carbons.items()}


def _named(mol, stand_in, forced):
    from .core import smiles_to_iupac

    smiles, position = _smiles_with_order(stand_in)
    token = FORCED_BRANCH_NAMES.set((stand_in.GetNumAtoms(), {position[i]: name for i, name in forced.items()}))
    try:
        return smiles_to_iupac(smiles)
    finally:
        FORCED_BRANCH_NAMES.reset(token)


def _bridge(mol, first, second):
    """The interior nitrogens of the plain N chain between acylated nitrogens `first` and `second`, else None."""
    interior = []
    previous, current = first, [n for n in mol.GetAtomWithIdx(first).GetNeighbors() if _plain_nitrogen(n)]
    if len(current) != 1:
        return None
    current = current[0].GetIdx()
    while current != second:
        atom = mol.GetAtomWithIdx(current)
        if atom.GetDegree() != 2 or atom.GetTotalNumHs() != 1:
            return None
        interior.append(current)
        following = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous]
        if not following or not _plain_nitrogen(mol.GetAtomWithIdx(following[0])):
            return None
        previous, current = current, following[0]
    return interior


def chain_amide_name(mol):
    """The name of a molecule with an acylated end of a nitrogen chain of three or more atoms, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return None
    ends = _acylated_ends(mol)
    if not ends:
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    if len(ends) == 2 and ends[0][2] != ends[1][0] and ends[1][2] != ends[0][0]:
        return _multiplicative(mol, graph, ends)
    if len(ends) != 1:
        return None
    nitrogen, root, partner = ends[0]
    try:
        donor = name_branch(graph, partner, nitrogen, halogens, aromatic, mol=mol)
    except UnsupportedStructure:
        return None
    stand_in, carbons = _stand_in(mol, [(nitrogen, partner)])
    try:
        name = _named(mol, stand_in, {carbons[partner]: donor})
    except UnsupportedStructure:
        return None
    return name if donor[0] in name and name.endswith("amide") else None


def _chalcogen_run(mol, start, behind):
    """The sulfur atoms of an unbranched run that starts at `start` and ends in SH, else None."""
    run, previous = [start], behind
    while True:
        atom = mol.GetAtomWithIdx(run[-1])
        if atom.GetAtomicNum() != 16 or atom.GetFormalCharge() or atom.IsInRing() or atom.GetIsotope():
            return None
        following = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous]
        if not following:
            return run if atom.GetTotalNumHs() == 1 else None
        if len(following) > 1 or mol.GetBondBetweenAtoms(run[-1], following[0]).GetBondTypeAsDouble() != 1.0:
            return None
        previous = run[-1]
        run.append(following[0])


def chain_ketone_name(mol):
    """P-58.3.2: a carbonyl group on a chain of three or more sulfur atoms is a ketone with a polysulfanyl prefix, not a
    pseudoester: 'phenyl(tetrasulfanyl)methanone'."""
    from ._hetero_prefixes import EXTENDED_PREFIXES

    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return None
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3:
            continue
        oxo = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
        sulfur = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 16]
        carbon = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(oxo) == 1 and len(sulfur) == 1 and len(carbon) == 1:
            run = _chalcogen_run(mol, sulfur[0].GetIdx(), atom.GetIdx())
            if run is not None and len(run) >= 3:
                found.append((atom.GetIdx(), run))
    if len(found) != 1:
        return None
    carbonyl, run = found[0]
    graph = adjacency(mol)
    token = EXTENDED_PREFIXES.set(True)
    try:
        donor = name_branch(graph, run[0], carbonyl, halogen_substituents(mol), mol=mol)
    except UnsupportedStructure:
        return None
    finally:
        EXTENDED_PREFIXES.reset(token)
    editable = Chem.RWMol(mol)
    for index in sorted(run, reverse=True):
        editable.RemoveAtom(index)
    carbonyl_index = carbonyl - sum(1 for r in run if r < carbonyl)
    ring = Chem.RWMol(Chem.CombineMols(editable.GetMol(), Chem.MolFromSmiles("c1ccccc1")))
    root = editable.GetNumAtoms()
    ring.AddBond(carbonyl_index, root, Chem.BondType.SINGLE)
    stand_in = ring.GetMol()
    try:
        Chem.SanitizeMol(stand_in)
        name = _named(mol, stand_in, {root: donor})
    except (UnsupportedStructure, Chem.rdchem.MolSanitizeException):
        return None
    return name if donor[0] in name and name.endswith("one") else None


def _multiplicative(mol, graph, ends):
    (n1, r1, _), (n2, r2, _) = ends
    interior = _bridge(mol, n1, n2)
    if not interior or len(interior) not in _LINKER:
        return None
    sides = [_side(graph, r1, n1), _side(graph, r2, n2)]
    if sides[0] & sides[1] or sides[0] | sides[1] | {n1, n2, *interior} != set(range(mol.GetNumAtoms())):
        return None
    keys = []
    for nitrogen, root, side in ((n1, r1, sides[0]), (n2, r2, sides[1])):
        part = Chem.RWMol(mol)
        for index in sorted(set(range(mol.GetNumAtoms())) - side, reverse=True):
            part.RemoveAtom(index)
        keys.append(Chem.MolToSmiles(part.GetMol()))
    if keys[0] != keys[1]:
        return None
    carbons = Chem.RWMol(mol)
    bridge = [n1, *interior, n2]
    for index in interior:
        carbons.GetAtomWithIdx(index).SetAtomicNum(6)
        carbons.GetAtomWithIdx(index).SetNoImplicit(False)
    stand_in = carbons.GetMol()
    try:
        Chem.SanitizeMol(stand_in)
        from .core import smiles_to_iupac

        name = smiles_to_iupac(Chem.MolToSmiles(stand_in))
    except (UnsupportedStructure, Chem.rdchem.MolSanitizeException):
        return None
    old, new = _LINKER[len(interior)]
    needle = f"({old})" if "-" in old else old
    if needle not in name or not name.endswith("amide"):
        return None
    return name.replace(needle, f"({new})" if "-" in new else new, 1)
