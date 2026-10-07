"""Naming of phosphonic acids (R-P(=O)(OH)2), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-67.1.1.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'phosphonic acid' is the preselected retained name for HP(O)(OH)2.
- P-67.1.1.2: when the hydrogen atom attached to the central phosphorus is
  substituted by an organyl group R, the substituent is cited directly as
  a prefix on the retained name -- 'ethylphosphonic acid' (PIN), *not*
  'ethanephosphonic acid'. The Blue Book explicitly rejects the
  sulfonic-acid-style suffix approach (which would give
  'benzenephosphonic acid') for exactly this reason: phosphonic acid's
  sibling phosphinic acid has *two* substitutable hydrogens on the central
  atom, and a suffix-on-parent-hydride approach would need extra letter
  locants to distinguish them, which the Blue Book calls unnecessarily
  cumbersome. Phosphonic acid itself has only one substitutable hydrogen,
  so there's never more than one R group and never a locant to omit or
  cite.
  Confirmed against real PubChem structures: 'methylphosphonic acid'
  (CID 13818, `CP(=O)(O)O`), 'phenylphosphonic acid' (CID 15295),
  'ethylphosphonic acid' (CID 204482), 'benzylphosphonic acid'
  (CID 81312) -- a ring substituent uses the exact same shape as an
  acyclic one, unlike the sulfonic/selenonic/telluronic acid modules
  elsewhere in this project, which need a separate plain-benzene-ring
  retained-name path.

Scope: a single phosphonic acid group (one phosphorus bonded to exactly
one non-oxygen substituent R, one double-bonded oxygen, and two hydroxyl
oxygens), with no other heteroatom anywhere in the molecule except the
group's own oxygens and, in R, any halogen. R is named via `name_branch`,
so it may be any substituent shape that helper already supports -- a
plain alkyl chain, a branched chain, a plain benzene ring, or a
halogenated variant of either. Explicitly out of scope (raise
`UnsupportedStructure`): a second phosphonic acid group, phosphinic acid
(two substituents on phosphorus -- see `_phosphinic_acid.py`),
coexistence with a more senior characteristic group, and any
functional-replacement/infix variant (phosphonous, phosphoric, etc.).
"""

from rdkit import Chem

from ._multiplicative_text import enclose
from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import name_branch

CENTER_STEMS = {15: "phosph", 33: "ars", 51: "stib"}


_SENIOR_ACIDS = [
    Chem.MolFromSmarts(smarts)
    for smarts in ("[CX3](=O)[OX2H1]", "[#16,#34,#52;X3,X4](=O)[OX2H1]", "[#34,#52;X4](=O)(=O)[OX2H1]")
]


def require_phosphorus_acid_scope(mol, phosphorus):
    """A single phosphorus whose acid outranks every other group of the molecule (P-41): carboxylic and the
    sulfur-group acids are senior to phosphonic and phosphinic acids, which become phosphono prefixes."""
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() in CENTER_STEMS and atom.GetIdx() != phosphorus.GetIdx():
            raise UnsupportedStructure("more than one phosphorus atom is not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR_ACIDS):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the phosphorus acid")


def _phosphonic_acid_phosphorus_atoms(mol):
    """Phosphorus atoms shaped like a phosphonic acid group: bonded to
    exactly one non-oxygen substituent, one double-bonded (terminal)
    oxygen, and two single-bonded hydroxyl oxygens (each terminal, one
    H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in CENTER_STEMS or atom.GetDegree() != 4 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(oxygens) != 3 or len(non_oxygens) != 1:
            continue
        double_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(hydroxyl_os) != 2:
            continue
        (double_o,) = double_os
        if double_o.GetDegree() != 1:
            continue
        if any(o.GetDegree() != 1 or o.GetTotalNumHs() != 1 for o in hydroxyl_os):
            continue
        matches.append(atom)
    return matches


def has_phosphonic_acid_shape(mol) -> bool:
    return bool(_phosphonic_acid_phosphorus_atoms(mol))


def name_phosphonic_acid(mol) -> str:
    phosphorus_atoms = _phosphonic_acid_phosphorus_atoms(mol)
    if len(phosphorus_atoms) != 1:
        raise UnsupportedStructure(
            "more than one phosphonic acid group is not supported yet"
        )
    (phosphorus,) = phosphorus_atoms

    require_phosphorus_acid_scope(mol, phosphorus)
    group_oxygens = {n.GetIdx() for n in phosphorus.GetNeighbors() if n.GetAtomicNum() == 8}

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    (root,) = [n for n in graph[phosphorus.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    name, is_compound = name_branch(graph, root, phosphorus.GetIdx(), halogens, aromatic_atoms, mol=mol)
    prefix = enclose(name) if is_compound else name
    return f"{prefix}{CENTER_STEMS[phosphorus.GetAtomicNum()]}onic acid"
