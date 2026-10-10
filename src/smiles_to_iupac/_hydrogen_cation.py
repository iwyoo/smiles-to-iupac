"""Acid salts drawn with explicit hydrons (P-65.6.2.3.1 method 2, P-65.6.2.3.2): the position of the acid hydrogen is not
given, so the hydrons are cited as the word 'hydrogen' (with 'di', 'tri') after the cations, and the anion is named fully
ionized ('sodium hydrogen 2-(carboxylatomethyl)benzoate').

The hydrons are stood in for by a metal cation that is absent from the salt, the salt is named, and that cation's word
is replaced by 'hydrogen' and moved behind the other cations."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._numerals import numerical_term

_STAND_INS = (("[Li+]", "lithium"), ("[Cs+]", "caesium"), ("[Rb+]", "rubidium"), ("[K+]", "potassium"), ("[Na+]", "sodium"))
_CATION_WORD = re.compile(r"^(?:di|tri|tetra)?(?:lithium|sodium|potassium|rubidium|caesium|magnesium|calcium|barium|strontium|azanium)$")


def hydrogen_salt_name(smiles, namer):
    parts = smiles.split(".")
    count = parts.count("[H+]")
    if not count or len(parts) < 2:
        return None
    others = [p for p in parts if p != "[H+]"]
    mol = Chem.MolFromSmiles(".".join(others))
    if mol is None or not any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()):
        return None
    present = {a.GetSymbol() for a in mol.GetAtoms() if a.GetFormalCharge() > 0}
    stand_in = next(((smi, word) for smi, word in _STAND_INS if Chem.MolFromSmiles(smi).GetAtomWithIdx(0).GetSymbol() not in present), None)
    if stand_in is None:
        return None
    smi, word = stand_in
    name = namer(".".join(others + [smi] * count))
    tokens = name.split(" ")
    cations = []
    for token in tokens:
        if _CATION_WORD.match(token):
            cations.append(token)
        else:
            break
    marker = (numerical_term(count) if count > 1 else "") + word
    if marker not in cations:
        raise UnsupportedStructure("the cation words of the salt name were not found")
    cations.remove(marker)
    hydrogen = (numerical_term(count) if count > 1 else "") + "hydrogen"
    return " ".join([*cations, hydrogen, *tokens[len(cations) + 1 :]])
