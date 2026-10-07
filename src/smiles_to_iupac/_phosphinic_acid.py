"""Naming of phosphinic acids (R2-P(=O)-OH), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-67.1.1.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'phosphinic acid' is the preselected retained name for H2P(O)(OH).
- P-67.1.1.2: like phosphonic acid (see `_phosphonic_acid.py`), when the
  hydrogen atoms attached to the central phosphorus are substituted by
  organyl groups, the substituents are cited directly as prefixes on the
  retained name -- PIN example: 'diethylphosphinic acid' (not
  'P-ethylethanephosphinic acid'). Unlike phosphonic acid, phosphinic
  acid has *two* substitutable hydrogens, so both substituents (possibly
  identical) are cited: `format_mononuclear_prefixes` already implements
  this exact "N substituents on one locant-free central atom" shape for
  `_phosphane.py`/`_phosphanone.py` (P-14.3.4.2(a) locants always
  omitted, P-16.5.1.3.1 alphabetization/parenthesization), so it's reused
  directly here rather than re-derived.
  Confirmed against real PubChem structures: 'dimethylphosphinic acid'
  (CID 76777, `CP(=O)(C)O`), 'diethylphosphinic acid' (CID 4186629),
  'diphenylphosphinic acid' (CID 15567), and the mixed-substituent case
  'ethyl(methyl)phosphinic acid' (CID 103893, `CCP(=O)(C)O`) -- no
  P-locants are used or needed (one P atom, nothing to number).

Scope: a single phosphinic acid group (one phosphorus bonded to exactly
two non-oxygen substituents R1/R2 -- possibly identical -- one
double-bonded oxygen, and one hydroxyl oxygen), with no other heteroatom
anywhere in the molecule except the group's own oxygens and, in R1/R2,
any halogen. Each R is named via `name_branch`, so it may be any
substituent shape that helper already supports (mirrors
`_phosphonic_acid.py`'s own scope exactly). Explicitly out of scope
(raise `UnsupportedStructure`): a second phosphinic acid group,
phosphonic acid (one substituent -- see `_phosphonic_acid.py`),
coexistence with a more senior characteristic group, and any
functional-replacement/infix variant.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._phosphonic_acid import CENTER_STEMS, require_phosphorus_acid_scope
from ._substituents import format_mononuclear_prefixes, name_branch



def _phosphinic_acid_phosphorus_atoms(mol):
    """Phosphorus atoms shaped like a phosphinic acid group: bonded to
    exactly two non-oxygen substituents, one double-bonded (terminal)
    oxygen, and one single-bonded hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in CENTER_STEMS or atom.GetDegree() != 4 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(oxygens) != 2 or len(non_oxygens) != 2:
            continue
        double_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(hydroxyl_os) != 1:
            continue
        (double_o,) = double_os
        if double_o.GetDegree() != 1:
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_phosphinic_acid_shape(mol) -> bool:
    return bool(_phosphinic_acid_phosphorus_atoms(mol))


def name_phosphinic_acid(mol) -> str:
    phosphorus_atoms = _phosphinic_acid_phosphorus_atoms(mol)
    if len(phosphorus_atoms) != 1:
        raise UnsupportedStructure(
            "more than one phosphinic acid group is not supported yet"
        )
    (phosphorus,) = phosphorus_atoms

    require_phosphorus_acid_scope(mol, phosphorus)
    group_oxygens = {n.GetIdx() for n in phosphorus.GetNeighbors() if n.GetAtomicNum() == 8}

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    roots = [n for n in graph[phosphorus.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    entries = [name_branch(graph, root, phosphorus.GetIdx(), halogens, aromatic_atoms, mol=mol) for root in roots]
    prefix = format_mononuclear_prefixes(entries)
    return f"{prefix}{CENTER_STEMS[phosphorus.GetAtomicNum()]}inic acid"
