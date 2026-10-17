"""Cyano groups on ring atoms, at least one on a ring heteroatom (P-66.5.1.1.3): the suffix 'carbonitrile', multiplied
for several groups, is used for -CN on a ring atom of any kind, 'piperidine-1-carbonitrile',
'piperazine-1,4-dicarbonitrile'. The groups are named as the carboxylic acids they rank beside, then the ending is
exchanged."""

from rdkit import Chem

from ._common import UnsupportedStructure

_NITRILE = Chem.MolFromSmarts("[NX1]#[CX2]")
_RING_N_CYANO = Chem.MolFromSmarts("[#7,#5,#13,#14,#15,#31,#32,#33,#49,#50,#51,#81,#82,#83;R;+0]-[CX2]#[NX1]")
_RING_CYANO = Chem.MolFromSmarts("[R;+0]-[CX2]#[NX1]")
_OTHER_SENIOR = Chem.MolFromSmarts("[#6,#16,#15](=[O,S,Se,Te,N])[O,N,S,F,Cl,Br,I]")


def has_ring_heteroatom_nitrile_shape(mol) -> bool:
    return bool(mol.HasSubstructMatch(_RING_N_CYANO)) and len(mol.GetSubstructMatches(_RING_CYANO)) == len(
        mol.GetSubstructMatches(_NITRILE)
    )


def _is_pseudoketone(mol, carbon_idx, double_idx):
    """A ring carbon whose heteroatom neighbours are all ring atoms joined by ring bonds: the carbonyl of a lactam, lactone
    or cyclic urea is a heterocyclic pseudoketone, junior to a nitrile (P-64.1.2.1, P-41)."""
    carbon = mol.GetAtomWithIdx(carbon_idx)
    if not carbon.IsInRing():
        return False
    for neighbor in carbon.GetNeighbors():
        if neighbor.GetIdx() == double_idx or neighbor.GetAtomicNum() == 6:
            continue
        if not neighbor.IsInRing() or not mol.GetBondBetweenAtoms(carbon_idx, neighbor.GetIdx()).IsInRing():
            return False
    return True


def _has_other_senior_group(mol) -> bool:
    return any(not _is_pseudoketone(mol, c, x) for c, x, _ in mol.GetSubstructMatches(_OTHER_SENIOR))


def name_ring_heteroatom_nitrile(mol) -> str:
    if not has_ring_heteroatom_nitrile_shape(mol) or _has_other_senior_group(mol):
        raise UnsupportedStructure("not a cyano group on a ring heteroatom alone")
    healed = Chem.RWMol(mol)
    matches = mol.GetSubstructMatches(_RING_CYANO)
    for _, carbon, nitrogen in matches:
        healed.GetAtomWithIdx(nitrogen).SetAtomicNum(8)
        healed.GetBondBetweenAtoms(carbon, nitrogen).SetBondType(Chem.BondType.DOUBLE)
        healed.AddBond(carbon, healed.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
    healed = healed.GetMol()
    Chem.SanitizeMol(healed)
    from .core import _smiles_to_iupac_unabridged

    name = _smiles_to_iupac_unabridged(Chem.MolToSmiles(healed))
    if not name.endswith("carboxylic acid") or name.count("carboxylic") != 1 or "carboxy" in name.replace("carboxylic", ""):
        raise UnsupportedStructure("the ring heteroatom acid analogue is not named by a carboxylic acid suffix")
    return name[: -len("carboxylic acid")] + "carbonitrile"
