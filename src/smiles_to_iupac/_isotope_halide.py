"""Isotopically modified acyl and sulfonyl halides (P-82.2.1, P-82.2.2, P-65.1.5, P-65.3): the halide takes the descriptor
of the acid it is derived from, because the nuclide goes before the same part of the name: '(2-13C)acetyl chloride'
from '(2-13C)acetic acid', '(13C)methanesulfonyl chloride' from '(13C)methanesulfonic acid'."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._numerals import multiplying_prefix

_GROUP = Chem.MolFromSmarts("[$([CX3](=O)),$([SX4](=O)=O)][F,Cl,Br,I]")
_OTHER_ACID_DERIVATIVE = Chem.MolFromSmarts("[CX3,SX4,PX4](=[O,S,N])[OX2,NX3,SX2]")
_HALIDE = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}


def has_labelled_halide_shape(mol) -> bool:
    return any(a.GetIsotope() for a in mol.GetAtoms()) and mol.HasSubstructMatch(_GROUP)


def name_labelled_halide(mol) -> str:
    groups = mol.GetSubstructMatches(_GROUP)
    halogens = {mol.GetAtomWithIdx(g[1]).GetAtomicNum() for g in groups}
    centres = {mol.GetAtomWithIdx(g[0]).GetAtomicNum() for g in groups}
    if len(halogens) != 1 or len(centres) != 1 or any(mol.GetAtomWithIdx(g[1]).GetIsotope() for g in groups):
        raise UnsupportedStructure("the acid halide groups of this molecule differ")
    if mol.HasSubstructMatch(_OTHER_ACID_DERIVATIVE):
        raise UnsupportedStructure("an acid, ester or amide group outranks the acid halide")
    editable = Chem.RWMol(mol)
    for _, halogen in groups:
        atom = editable.GetAtomWithIdx(halogen)
        atom.SetAtomicNum(8)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(1)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    from .core import smiles_to_iupac

    name = smiles_to_iupac(Chem.MolToSmiles(acid))
    if name.endswith("carboxylic acid"):
        stem = name[: -len("ic acid")][: -len("carboxyl")] + "carbonyl"
    elif name.endswith("ic acid") and name.count(" acid") == 1:
        stem = name[: -len("ic acid")] + "yl"
    else:
        raise UnsupportedStructure("the acid analogue of this halide has no acid ending to exchange")
    if _every_hydrogen_is_one_nuclide(mol):
        stem = re.sub(r"^\([\d,]+-(\d+H\d+)\)", r"(\1)", stem)
    count = len(groups)
    return f"{stem} {multiplying_prefix(count) if count > 1 else ''}{_HALIDE[next(iter(halogens))]}"


def _every_hydrogen_is_one_nuclide(mol):
    """The locants of a hydrogen nuclide are dropped when it replaces all the hydrogens of the parent (P-82.2.1)."""
    explicit = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    return bool(explicit) and len({a.GetIsotope() for a in explicit}) == 1 and not any(a.GetTotalNumHs() for a in mol.GetAtoms())
