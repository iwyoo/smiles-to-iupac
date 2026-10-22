"""Imidoyl and hydrazonoyl cations (P-73.2.3.1, P-65.1.7): the carboximidic or carbohydrazonic acid name with
'ic acid' replaced by 'oylium'."""

from rdkit import Chem

from ._common import UnsupportedStructure

_ENDINGS = (("imidic acid", "imidoylium"), ("hydrazonic acid", "hydrazonoylium"))


def _center(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    atom = charged[0]
    if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 1 or atom.GetIsotope() or atom.IsInRing():
        return None
    neighbors = list(atom.GetNeighbors())
    double = [
        n for n in neighbors if n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if len(double) != 1 or atom.GetNumRadicalElectrons() or len(neighbors) - 1 + atom.GetTotalNumHs() != 1:
        return None
    if any(n.GetAtomicNum() != 6 for n in neighbors if n.GetIdx() != double[0].GetIdx()):
        return None
    return atom


def has_imidoyl_cation_shape(mol) -> bool:
    return _center(mol) is not None


def name_imidoyl_cation(mol) -> str:
    from .core import smiles_to_iupac

    center = _center(mol)
    if center is None:
        raise UnsupportedStructure("not an imidoyl cation")
    editable = Chem.RWMol(mol)
    atom = editable.GetAtomWithIdx(center.GetIdx())
    atom.SetFormalCharge(0)
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(0)
    editable.AddBond(center.GetIdx(), editable.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    name = smiles_to_iupac(Chem.MolToSmiles(acid))
    for tail, replacement in _ENDINGS:
        if name.endswith(tail):
            return name[: -len(tail)] + replacement
    raise UnsupportedStructure("this imidoyl cation has no imidic or hydrazonic acid name")
