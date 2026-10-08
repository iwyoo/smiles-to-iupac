"""A cyano group on a ring nitrogen (P-66.5.1.1.3): the suffix 'carbonitrile' is used for -CN on a ring atom of any
kind, 'piperidine-1-carbonitrile'. The group is named as the carboxylic acid it ranks beside, then the ending is
exchanged."""

from rdkit import Chem

from ._common import UnsupportedStructure

_NITRILE = Chem.MolFromSmarts("[NX1]#[CX2]")
_RING_N_CYANO = Chem.MolFromSmarts("[#7;R;X3;+0]-[CX2]#[NX1]")
_OTHER_SENIOR = Chem.MolFromSmarts("[#6,#16,#15](=[O,S,Se,Te,N])[O,N,S,F,Cl,Br,I]")


def has_ring_heteroatom_nitrile_shape(mol) -> bool:
    return len(mol.GetSubstructMatches(_RING_N_CYANO)) == 1 and len(mol.GetSubstructMatches(_NITRILE)) == 1


def name_ring_heteroatom_nitrile(mol) -> str:
    if not has_ring_heteroatom_nitrile_shape(mol) or mol.HasSubstructMatch(_OTHER_SENIOR):
        raise UnsupportedStructure("not a cyano group on a ring nitrogen alone")
    (_, carbon, nitrogen) = mol.GetSubstructMatch(_RING_N_CYANO)
    healed = Chem.RWMol(mol)
    healed.GetAtomWithIdx(nitrogen).SetAtomicNum(8)
    healed.GetBondBetweenAtoms(carbon, nitrogen).SetBondType(Chem.BondType.DOUBLE)
    hydroxyl = healed.AddAtom(Chem.Atom(8))
    healed.AddBond(carbon, hydroxyl, Chem.BondType.SINGLE)
    healed = healed.GetMol()
    Chem.SanitizeMol(healed)
    from .core import _smiles_to_iupac_unabridged

    name = _smiles_to_iupac_unabridged(Chem.MolToSmiles(healed))
    if not name.endswith("-carboxylic acid") or name.count("carboxylic") != 1:
        raise UnsupportedStructure("the ring nitrogen acid analogue is not named by a carboxylic acid suffix")
    return name[: -len("carboxylic acid")] + "carbonitrile"
