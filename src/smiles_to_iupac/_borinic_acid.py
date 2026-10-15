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
group's own oxygen and, in R1/R2, nitrogen, sulfur or any halogen. Each R is named via
`name_branch`, so it may be any substituent shape that helper already
supports (mirrors `_boronic_acid.py`/`_phosphinic_acid.py`'s own scope
exactly). Explicitly out of scope (raise `UnsupportedStructure`): a
second borinic acid group, boronic acid (one substituent -- see
`_boronic_acid.py`), coexistence with a more senior characteristic
group, and any chalcogen-replacement analogue.
"""

from rdkit import Chem

from ._acid_groups import acid_group_at
from ._hetero_prefixes import ACIDS_SENIOR_TO_BORON, boron_acid_parent
from ._parent_choice import senior_parent_name
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import cited_stereo_around, format_mononuclear_prefixes, name_branch

_BORON = 5
BORON_ACID_ELEMENTS = (1, 6, 7, 8, 16, _BORON, *HALOGEN_PREFIXES)


def _borinic_acid_boron_atoms(mol):
    """Boron atoms shaped like a borinic acid group: bonded to exactly two
    non-oxygen substituents and one single-bonded hydroxyl oxygen
    (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _BORON or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0 or atom.IsInRing():
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
    from ._boronic_acid import _boronic_acid_boron_atoms

    boron_atoms = _borinic_acid_boron_atoms(mol)
    if not boron_atoms:
        raise UnsupportedStructure("no borinic acid group")
    if _boronic_acid_boron_atoms(mol):
        raise UnsupportedStructure("a boronic acid group outranks the borinic acid (P-42.3)")
    borons = {a.GetIdx() for a in boron_atoms}

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _BORON and atom.GetIdx() not in borons:
            raise UnsupportedStructure("more than one boron atom is not supported yet")
        if atomic_num not in BORON_ACID_ELEMENTS:
            raise UnsupportedStructure(
                "heteroatoms other than the borinic acid's own boron/"
                "oxygen, nitrogen, sulfur and a halogen substituent are not supported yet"
            )

    if any(mol.HasSubstructMatch(query) for query in ACIDS_SENIOR_TO_BORON) or any(
        acid_group_at(mol, a.GetIdx()) is not None for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the borinic acid")

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    if len(borons) == 1:
        return _name_at(mol, boron_atoms[0])[1]
    return senior_parent_name(mol, sorted(borons), lambda idx: _name_at(mol, mol.GetAtomWithIdx(idx)))


def _name_at(mol, boron):
    """(number of substituent prefixes, name) with the borinic acid group at `boron` as the parent."""
    group_oxygens = {n.GetIdx() for n in boron.GetNeighbors() if n.GetAtomicNum() == 8}
    graph = adjacency(mol)
    roots = [n for n in graph[boron.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    with boron_acid_parent(), cited_stereo_around(mol, boron.GetIdx()):
        entries = [name_branch(graph, root, boron.GetIdx(), halogens, aromatic_atoms, mol=mol) for root in roots]
    return len(entries), f"{format_mononuclear_prefixes(entries)}borinic acid"
