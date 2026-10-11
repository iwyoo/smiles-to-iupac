"""Anions of the organyl esters of thio-, dithio- and trithiocarbonic acids (P-65.2.1.2, P-65.6.3.3.7.2.1, P-72.2.2.2.1): the
anion of 'O-ethyl hydrogen carbonodithioate' is 'O-ethyl carbonodithioate', the word 'hydrogen' standing for the proton
that the anionic chalcogen has lost. The partly esterified acid is named first, so the O/S locants of the organyl
groups follow its rules."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure

_CARBONIC_ANION = Chem.MolFromSmarts("[CX3](=[O,S,Se,Te])([O,S,Se,Te])[O,S,Se,Te;X1;-1]")


def has_carbonic_ester_anion_shape(mol) -> bool:
    if len(Chem.GetMolFrags(mol)) != 1 or sum(a.GetFormalCharge() for a in mol.GetAtoms()) != -1:
        return False
    return sum(1 for a in mol.GetAtoms() if a.GetFormalCharge()) == 1 and mol.HasSubstructMatch(_CARBONIC_ANION)


def name_carbonic_ester_anion(mol) -> str:
    if not has_carbonic_ester_anion_shape(mol):
        raise UnsupportedStructure("not the anion of a carbonic acid ester")
    (_, _, _, anionic), *_ = mol.GetSubstructMatches(_CARBONIC_ANION)
    protonated = Chem.RWMol(mol)
    atom = protonated.GetAtomWithIdx(anionic)
    atom.SetFormalCharge(0)
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(1)
    acid = protonated.GetMol()
    Chem.SanitizeMol(acid)
    from .core import smiles_to_iupac

    name = smiles_to_iupac(Chem.MolToSmiles(acid))
    if name.startswith(("hydrogen ", "dihydrogen ")) or " hydrogen " not in name:
        raise UnsupportedStructure("the acid of this anion is not named as an ester with a hydrogen word")
    return re.sub(r" hydrogen ", " ", name, count=1)
