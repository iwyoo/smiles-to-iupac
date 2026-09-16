"""Naming of boronic acids (R-B(OH)2), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-68.1.4.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'boronic acid' is the preselected retained name for HB(OH)2, and when
  the hydrogen atom attached to boron is substituted by an organyl group
  R, the substituent is cited directly as a prefix on the retained name --
  'methylboronic acid' for CH3-B(OH)2 (PIN), not the substitutive
  'methaneboronic acid' (see P-67.1.1.2, the same rule `_phosphonic_acid.py`
  already applies for the phosphorus analogue -- identical shape, boron in
  place of phosphorus, no P=O).
  Confirmed against real PubChem structures: 'methylboronic acid'
  (CID 139377, `CB(O)O`), 'phenylboronic acid' (CID 66827),
  'butylboronic acid' (CID 20479), 'benzylboronic acid' (CID 11320956) --
  a ring substituent uses the exact same shape as an acyclic one.

Scope: a single boronic acid group (one boron bonded to exactly one
non-oxygen substituent R and two hydroxyl oxygens -- boron has no
double-bonded oxygen here, unlike phosphonic acid's P=O), with no other
heteroatom anywhere in the molecule except the group's own oxygens and,
in R, any halogen. R is named via `name_branch`, so it may be any
substituent shape that helper already supports -- a plain alkyl chain, a
branched chain, a plain benzene ring, or a halogenated variant of either
(mirrors `_phosphonic_acid.py`'s own scope exactly). Explicitly out of
scope (raise `UnsupportedStructure`): a second boronic acid group, borinic
acid (two substituents on boron -- see `_borinic_acid.py`), coexistence
with a more senior characteristic group, and any chalcogen-replacement
analogue (thioboronic acid etc.).
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import name_branch

_BORON = 5


def _boronic_acid_boron_atoms(mol):
    """Boron atoms shaped like a boronic acid group: bonded to exactly one
    non-oxygen substituent and two single-bonded hydroxyl oxygens (each
    terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _BORON or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(oxygens) != 2 or len(non_oxygens) != 1:
            continue
        if any(
            mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() != 1.0 for o in oxygens
        ):
            continue
        if any(o.GetDegree() != 1 or o.GetTotalNumHs() != 1 for o in oxygens):
            continue
        matches.append(atom)
    return matches


def has_boronic_acid_shape(mol) -> bool:
    return bool(_boronic_acid_boron_atoms(mol))


def name_boronic_acid(mol) -> str:
    boron_atoms = _boronic_acid_boron_atoms(mol)
    if len(boron_atoms) != 1:
        raise UnsupportedStructure("more than one boronic acid group is not supported yet")
    (boron,) = boron_atoms

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _BORON and atom.GetIdx() != boron.GetIdx():
            raise UnsupportedStructure("more than one boron atom is not supported yet")
        if atomic_num not in (1, 6, 8, _BORON, *HALOGEN_PREFIXES):
            raise UnsupportedStructure(
                "heteroatoms other than the boronic acid's own boron/"
                "oxygens and a halogen substituent are not supported yet"
            )

    group_oxygens = {n.GetIdx() for n in boron.GetNeighbors() if n.GetAtomicNum() == 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetIdx() not in group_oxygens:
            raise UnsupportedStructure(
                "an oxygen atom not part of the boronic acid's own "
                "B(OH)2 group is out of scope for this module"
            )

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    (root,) = [n for n in graph[boron.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    name, is_compound = name_branch(graph, root, boron.GetIdx(), halogens, aromatic_atoms, mol=mol)
    prefix = f"({name})" if is_compound else name
    return f"{prefix}boronic acid"
