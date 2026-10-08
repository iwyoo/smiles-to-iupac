"""Linear peptides named by the acyl groups of their N-terminal residues and the C-terminal amino acid (P-103.3.2, P-103.3.4).

The molecule is cut at every eupeptide bond; each piece is named as the free amino acid and the names of all but the
C-terminal one change their ending to 'yl' (P-103.2.5).
"""

import re
from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure

_PEPTIDE_BOND = Chem.MolFromSmarts("[CX3](=O)(-[CX4]-[NX3])-[NX3;!$(N=*)]-[CX4]-[CX3]=O")
_SPECIAL_ACYL = {
    "cysteine": "cysteinyl",
    "asparagine": "asparaginyl",
    "glutamine": "glutaminyl",
    "tryptophan": "tryptophyl",
    "aspartic acid": "aspartyl",
    "glutamic acid": "glutamyl",
}
_RESIDUES = {
    "alanine", "arginine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine",
    "proline", "serine", "threonine", "valine", "tyrosine", "glycine", *_SPECIAL_ACYL,
}
_RESIDUE_NAME = re.compile(r"^(?:([LD])-)?(allo)?([a-z ]+)$")


def _acyl(base):
    return _SPECIAL_ACYL.get(base) or base[: -len("ine")] + "yl"


def _cut(mol):
    cuts = {(m[0], m[4]) for m in mol.GetSubstructMatches(_PEPTIDE_BOND)}
    editable = Chem.RWMol(mol)
    for carbon, nitrogen in cuts:
        editable.RemoveBond(carbon, nitrogen)
        editable.AddBond(carbon, editable.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    pieces = editable.GetMol()
    Chem.SanitizeMol(pieces)
    return cuts, pieces


def _sequence(cuts, fragments):
    """Fragment indices from the N-terminal residue to the C-terminal one, None unless they form one chain."""
    owner = {a: i for i, atoms in enumerate(fragments) for a in atoms}
    following = {owner[c]: owner[n] for c, n in cuts}
    if len(following) != len(cuts) or len(set(following.values())) != len(cuts):
        return None
    first = [i for i in range(len(fragments)) if i not in following.values()]
    if len(first) != 1:
        return None
    order = [first[0]]
    while order[-1] in following:
        order.append(following[order[-1]])
    return order if len(order) == len(fragments) else None


@lru_cache(maxsize=256)
def _name(smiles):
    from .core import smiles_to_iupac

    mol = Chem.MolFromSmiles(smiles)
    cuts, pieces = _cut(mol)
    fragments = Chem.GetMolFrags(pieces)
    if not cuts or len(fragments) != len(cuts) + 1:
        return None
    order = _sequence(cuts, fragments)
    if order is None:
        return None
    molecules = Chem.GetMolFrags(pieces, asMols=True)
    words = []
    for position, index in enumerate(order):
        try:
            residue = smiles_to_iupac(Chem.MolToSmiles(molecules[index]))
        except (UnsupportedStructure, ValueError):
            return None
        found = _RESIDUE_NAME.match(residue)
        if not found or found.group(3) not in _RESIDUES:
            return None
        descriptor, allo, base = found.groups()
        if descriptor is None and base != "glycine":
            descriptor = "ξ"
        word = (allo or "") + (_acyl(base) if position < len(order) - 1 else base)
        words.append((f"{descriptor}-" if descriptor else "") + word)
    return words[0] + "".join(("-" if w[1:2] == "-" else "") + w for w in words[1:])


def has_peptide_shape(mol) -> bool:
    return mol.HasSubstructMatch(_PEPTIDE_BOND) and _name(Chem.MolToSmiles(mol)) is not None


def name_peptide(mol) -> str:
    return _name(Chem.MolToSmiles(mol))
