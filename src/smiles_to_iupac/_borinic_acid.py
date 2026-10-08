"""Naming of borinic acids (R2-B-OH), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-68.1.4.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'borinic acid' is the preselected retained name for H2B(OH). Like
  phosphinic acid (see `_phosphinic_acid.py`, the phosphorus analogue --
  same "two substitutable hydrogens on one locant-free central atom"
  shape), both substituents are cited directly as prefixes on the
  retained name, using `format_mononuclear_prefixes` for the
  alphabetization/parenthesization/multiplying-prefix logic
  (P-16.5.1.3.1).
  Confirmed against real PubChem structures: 'dimethylborinic acid'
  (CID 5326208, `B(C)(C)O`), 'diethylborinic acid' (CID 545119),
  'diphenylborinic acid' (CID 17498), and the mixed-substituent case
  'ethyl(methyl)borinic acid' (CID 22096092, `B(C)(CC)O`) -- no
  B-locants are used or needed (one B atom, nothing to number).

Scope: a single borinic acid group (one boron bonded to exactly two
non-oxygen substituents R1/R2 -- possibly identical -- and one hydroxyl
oxygen; boron has no double-bonded oxygen here, unlike phosphinic acid's
P=O), with no other heteroatom anywhere in the molecule except the
group's own oxygen and, in R1/R2, any halogen. Each R is named via
`name_branch`, so it may be any substituent shape that helper already
supports (mirrors `_boronic_acid.py`/`_phosphinic_acid.py`'s own scope
exactly). Explicitly out of scope (raise `UnsupportedStructure`): a
second borinic acid group, boronic acid (one substituent -- see
`_boronic_acid.py`), coexistence with a more senior characteristic
group, and any chalcogen-replacement analogue.
"""

from rdkit import Chem

from ._hetero_prefixes import ACIDS_SENIOR_TO_BORON
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_mononuclear_prefixes, name_branch

_BORON = 5


def _borinic_acid_boron_atoms(mol):
    """Boron atoms shaped like a borinic acid group: bonded to exactly two
    non-oxygen substituents and one single-bonded hydroxyl oxygen
    (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _BORON or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(oxygens) != 1 or len(non_oxygens) != 2:
            continue
        (hydroxyl_o,) = oxygens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), hydroxyl_o.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_borinic_acid_shape(mol) -> bool:
    return bool(_borinic_acid_boron_atoms(mol))


def name_borinic_acid(mol) -> str:
    boron_atoms = _borinic_acid_boron_atoms(mol)
    if len(boron_atoms) != 1:
        raise UnsupportedStructure("more than one borinic acid group is not supported yet")
    (boron,) = boron_atoms

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _BORON and atom.GetIdx() != boron.GetIdx():
            raise UnsupportedStructure("more than one boron atom is not supported yet")
        if atomic_num not in (1, 6, 8, _BORON, *HALOGEN_PREFIXES):
            raise UnsupportedStructure(
                "heteroatoms other than the borinic acid's own boron/"
                "oxygen and a halogen substituent are not supported yet"
            )

    group_oxygens = {n.GetIdx() for n in boron.GetNeighbors() if n.GetAtomicNum() == 8}
    if any(mol.HasSubstructMatch(query) for query in ACIDS_SENIOR_TO_BORON):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the borinic acid")

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    roots = [n for n in graph[boron.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    entries = [name_branch(graph, root, boron.GetIdx(), halogens, aromatic_atoms, mol=mol) for root in roots]
    prefix = format_mononuclear_prefixes(entries)
    return f"{prefix}borinic acid"
