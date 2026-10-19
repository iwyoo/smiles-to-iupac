"""A nuclide on the O-H or N-H of a boron, phosphorus, arsenic, antimony, sulfur, carbonic, oxalic, cyanic or carbamic acid
group or of cyanamide (P-82.2.1, P-82.6.2): the compound is named without it and the descriptor, with the letter locant of
the modified atom, is inserted before the parent, after its prefixes: 'methyl(O-2H)phosphonic acid'; the hydrogen of a partial
ester is the word 'hydrogen' that carries it, the modified one cited first (P-82.2.2.1): 'methyl (2H)hydrogen hydrogen phosphate'."""

import re
from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure
from ._isotope_labels import restore_labels, split_isotopes

_ACID_PARENT = re.compile(
    r"(?:(?:phosph|ars|stib)[a-z]*(?:ic|ous)|sulfur\w*(?:ic|ous)|sulfamic|bor(?:on|in)?ic|carbamic|carbonic|oxalic|cyanic) acid$|cyanamide$"
)
_ACID_HYDROGEN_WORD = re.compile(r"(?<![a-z])(di)?hydrogen(?= )")


def _is_oxoacid_centre(atom):
    """Whether a hydroxy or amino group on `atom` belongs to an oxoacid: a terminal =O on carbon, phosphorus, sulfur,
    arsenic or antimony, a boron atom, or a carbon triply bonded to nitrogen."""
    z = atom.GetAtomicNum()
    if z == 5:
        return True
    if z not in (6, 15, 16, 33, 51):
        return False
    return any(
        (b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 8 and b.GetOtherAtom(atom).GetDegree() == 1)
        or (z == 6 and b.GetBondTypeAsDouble() == 3.0 and b.GetOtherAtom(atom).GetAtomicNum() == 7)
        for b in atom.GetBonds()
    )


def _is_carbamic_or_cyanamide_centre(carbon):
    """A carbon that bears the N-H of carbamic acid (a hydroxy group beside =O) or of cyanamide (C#N); an amide is no acid."""
    if carbon.GetAtomicNum() != 6:
        return False
    return any(
        (b.GetBondTypeAsDouble() == 3.0 and b.GetOtherAtom(carbon).GetAtomicNum() == 7)
        or (b.GetBondTypeAsDouble() == 1.0 and b.GetOtherAtom(carbon).GetAtomicNum() == 8 and b.GetOtherAtom(carbon).GetTotalNumHs() > 0)
        for b in carbon.GetBonds()
    )


def _labelled_group_atom(mol):
    """(mol without the nuclide of the acid hydrogen, its atom, the label entry, the atoms with the same label) for a
    structure whose only acid hydrogens modified are one nuclide on oxoacid O-H groups (or one N-H); the other nuclides
    stay where they are."""
    split = split_isotopes(mol)
    if split is None or len(Chem.GetMolFrags(mol)) != 1:
        return None
    clean, labels, _ = split
    acid = {}
    for index, entry in labels.items():
        atom = clean.GetAtomWithIdx(index)
        if entry["skeleton"] or sum(entry["H"].values()) != 1 or atom.GetAtomicNum() not in (7, 8):
            continue
        centres = [n for n in atom.GetNeighbors() if _is_oxoacid_centre(n)]
        if len(centres) == 1 and (atom.GetAtomicNum() == 8 or _is_carbamic_or_cyanamide_centre(centres[0])):
            acid[index] = entry
    if not acid or len({tuple(entry["H"]) for entry in acid.values()}) != 1:
        return None
    if len(acid) > 1 and any(clean.GetAtomWithIdx(i).GetAtomicNum() != 8 for i in acid):
        return None
    first = min(acid)
    editable = Chem.RWMol(clean)
    restore_labels(editable, clean, labels, skip=set(acid))
    return editable.GetMol(), clean.GetAtomWithIdx(first), acid[first], sorted(acid)


@lru_cache(maxsize=256)
def _unlabelled_name(smiles):
    from .core import smiles_to_iupac

    try:
        return smiles_to_iupac(smiles)
    except UnsupportedStructure:
        return None


def _named_group(mol):
    """The `_labelled_group_atom` result with the name of the structure without the acid hydrogen label, or None when
    that name holds neither an acid parent nor the word 'hydrogen' that can carry the label."""
    found = _labelled_group_atom(mol)
    if found is None:
        return None
    name = _unlabelled_name(Chem.MolToSmiles(found[0]))
    if name is None or not (_ACID_PARENT.search(name) or _ACID_HYDROGEN_WORD.search(name)):
        return None
    return found, name


def has_isotope_acid_hydrogen_shape(mol) -> bool:
    return _named_group(mol) is not None


def name_isotope_acid_hydrogen(mol) -> str:
    named = _named_group(mol)
    if named is None:
        raise UnsupportedStructure("no labelled hydrogen on an oxoacid group")
    (_, atom, entry, modified), name = named
    (nuclide,) = entry["H"]
    match = _ACID_PARENT.search(name)
    if match is None:
        return _hydrogen_word_citation(name, nuclide, len(modified), atom)
    if len(modified) > 1 or re.search(r"(?<![a-z])(?:di|tri|tetra)$", name[: match.start()]):
        raise UnsupportedStructure("several modified acid hydrogens are not cited on this acid yet")
    count = "1" if atom.GetTotalNumHs() > 1 else ""
    return f"{name[: match.start()]}({atom.GetSymbol()}-{nuclide}{count}){name[match.start():]}"


def _hydrogen_word_citation(name, nuclide, modified, atom):
    """The ester name with the nuclide in front of the word 'hydrogen' that holds the modified atom; a nuclide cited
    separately from the unmodified hydrogen, which follows it (P-82.2.2.1, P-82.2.4)."""
    word = _ACID_HYDROGEN_WORD.search(name)
    present = 2 if word.group(1) else 1
    if modified > present:
        raise UnsupportedStructure("more modified acid hydrogens than the acid hydrogens of the ester")
    if modified == present:
        citation = f"({nuclide}{modified if modified > 1 else ''})" + ("dihydrogen" if present == 2 else "hydrogen")
    else:
        citation = f"({nuclide})hydrogen hydrogen"
    return name[: word.start()] + citation + name[word.end():]
