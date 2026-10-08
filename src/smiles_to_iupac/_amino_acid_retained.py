"""Retained names of the less common amino acids of Table 10.5 whose skeleton is not an alpha-amino acid with a side chain (P-103.1.1.2).

The stereocentres (atom-mapped in the templates) give the D or L descriptor of P-103.1.3.1: S is L, except where a sulfur
on C-3 outranks the carboxy group, which reverses the correspondence.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

_TEMPLATES = (
    ("β-alanine", "NCCC(=O)O", False),
    ("sarcosine", "CNCC(=O)O", False),
    ("5-oxoproline", "O=C1CC[CH:1](N1)C(=O)O", False),
    ("homoserine lactone", "N[CH:1]1CCOC1=O", False),
    ("thyroxine", "N[CH:1](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O", False),
    ("cystine", "N[CH:1](CSSC[CH:2](N)C(=O)O)C(=O)O", True),
    ("lanthionine", "N[CH:1](CSC[CH:2](N)C(=O)O)C(=O)O", True),
)
_COMPILED = tuple((name, Chem.MolFromSmiles(smiles), sulfur) for name, smiles, sulfur in _TEMPLATES)


def _match(mol):
    for name, template, sulfur in _COMPILED:
        if any(atom.GetFormalCharge() or atom.GetIsotope() for atom in mol.GetAtoms()):
            return None
        if mol.GetNumAtoms() == template.GetNumAtoms() and mol.GetNumBonds() == template.GetNumBonds():
            found = mol.GetSubstructMatch(template)
            if found:
                centers = [found[a.GetIdx()] for a in template.GetAtoms() if a.GetAtomMapNum()]
                return name, centers, sulfur
    return None


def _descriptor(mol, centers, sulfur):
    """'L', 'D', '' when no centre is specified, None when they disagree or only some are specified."""
    if not centers:
        return ""
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    labels = [probe.GetAtomWithIdx(c).GetPropsAsDict().get("_CIPCode") for c in centers]
    if all(label is None for label in labels):
        return ""
    if None in labels or len(set(labels)) != 1:
        return None
    mapping = {"R": "L", "S": "D"} if sulfur else {"S": "L", "R": "D"}
    return mapping[labels[0]]


def has_retained_amino_acid_shape(mol) -> bool:
    found = _match(mol)
    return found is not None and _descriptor(mol, found[1], found[2]) is not None


def name_retained_amino_acid(mol) -> str:
    name, centers, sulfur = _match(mol)
    descriptor = _descriptor(mol, centers, sulfur)
    return f"{descriptor}-{name}" if descriptor else name
