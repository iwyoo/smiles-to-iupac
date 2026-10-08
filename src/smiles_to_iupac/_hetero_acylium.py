"""Acylium cations of sulfinic, sulfonic and phosphorus oxoacids (P-73.2.3.1): the acid name with 'ic acid' replaced by
'ylium' (ethenesulfinylium, benzenesulfonylium), 'inic acid' by 'inoylium' (dimethylphosphinoylium) and the
phosphonic acid by 'phosphonobis(ylium)'."""

from rdkit import Chem

from ._common import UnsupportedStructure

_ENDINGS = (
    ("sulfinic acid", "sulfinylium"),
    ("sulfonic acid", "sulfonylium"),
    ("phosphinic acid", "phosphinoylium"),
    ("phosphonic acid", "phosphonobis(ylium)"),
)


def _center(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or len(Chem.GetMolFrags(mol)) != 1 or charged[0].GetIsotope():
        return None
    atom = charged[0]
    z, charge = atom.GetAtomicNum(), atom.GetFormalCharge()
    oxo = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
    oxo_indices = {n.GetIdx() for n in oxo}
    others = [n for n in atom.GetNeighbors() if n.GetIdx() not in oxo_indices]
    if atom.IsInRing() or atom.GetTotalNumHs() or any(n.GetAtomicNum() != 6 for n in others):
        return None
    shape = (z, charge, len(oxo), len(others))
    if shape in ((16, 1, 1, 1), (16, 1, 2, 1), (15, 1, 1, 2), (15, 2, 1, 1)):
        if all(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0 for n in oxo):
            return atom
    return None


def has_hetero_acylium_shape(mol) -> bool:
    return _center(mol) is not None


def name_hetero_acylium(mol) -> str:
    from .core import smiles_to_iupac

    center = _center(mol)
    if center is None:
        raise UnsupportedStructure("not an acylium cation of a sulfur or phosphorus oxoacid")
    editable = Chem.RWMol(mol)
    atom = editable.GetAtomWithIdx(center.GetIdx())
    missing = atom.GetFormalCharge()
    atom.SetFormalCharge(0)
    atom.SetNoImplicit(False)
    for _ in range(missing):
        hydroxy = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(center.GetIdx(), hydroxy, Chem.BondType.SINGLE)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    name = smiles_to_iupac(Chem.MolToSmiles(acid))
    for tail, replacement in _ENDINGS:
        if name.endswith(tail):
            return name[: -len(tail)] + replacement
    raise UnsupportedStructure("this acylium has no oxoacid name to convert")
