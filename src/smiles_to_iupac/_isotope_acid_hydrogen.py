"""A nuclide on the O-H or N-H of a phosphorus, arsenic, antimony, sulfur or carbamic acid group (P-82.2.1, P-82.6.2): the
acid is named without it and the descriptor, with the letter locant of the modified atom, is inserted before the acid
parent, after its prefixes: 'methyl(O-2H)phosphonic acid'."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._isotope_labels import split_isotopes

_ACID_PARENT = re.compile(
    r"(?:(?:phosph|ars|stib)(?:or|on|in)(?:ic|ous)|sulfur(?:ic|ous)|carbamic) acid$"
)


def _labelled_group_atom(mol):
    split = split_isotopes(mol)
    if split is None or len(Chem.GetMolFrags(mol)) != 1:
        return None
    clean, labels, _ = split
    if len(labels) != 1:
        return None
    ((index, entry),) = labels.items()
    atom = clean.GetAtomWithIdx(index)
    if entry["skeleton"] or sum(entry["H"].values()) != 1 or atom.GetAtomicNum() not in (7, 8):
        return None
    centres = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() in (6, 15, 16, 33, 51)
        and any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(n).GetAtomicNum() == 8 and b.GetOtherAtom(n).GetDegree() == 1
            for b in n.GetBonds()
        )
    ]
    if len(centres) != 1 or (atom.GetAtomicNum() == 7 and centres[0].GetAtomicNum() != 6):
        return None
    return clean, atom, entry


def has_isotope_acid_hydrogen_shape(mol) -> bool:
    return _labelled_group_atom(mol) is not None


def name_isotope_acid_hydrogen(mol) -> str:
    from .core import smiles_to_iupac

    found = _labelled_group_atom(mol)
    if found is None:
        raise UnsupportedStructure("no labelled hydrogen on an oxoacid group")
    clean, atom, entry = found
    name = smiles_to_iupac(Chem.MolToSmiles(clean))
    match = _ACID_PARENT.search(name)
    if match is None:
        raise UnsupportedStructure("the parent acid of this isotopically modified structure is not delimited")
    (nuclide,) = entry["H"]
    count = "1" if atom.GetTotalNumHs() > 1 else ""
    return f"{name[: match.start()]}({atom.GetSymbol()}-{nuclide}{count}){name[match.start():]}"
