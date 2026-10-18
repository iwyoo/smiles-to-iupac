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

import contextlib

from rdkit import Chem

from ._hetero_prefixes import EXTENDED_PREFIXES
from ._multiplicative_text import enclose
from ._common import UnsupportedStructure, adjacency, halogen_substituents, nitro_group_atoms
from ._substituents import centre_nuclides, cited_stereo_around, name_branch

CENTER_STEMS = {15: "phosph", 33: "ars", 51: "stib"}


_SENIOR_ACIDS = [
    Chem.MolFromSmarts(smarts)
    for smarts in ("[CX3](=O)[OX2H1]", "[#16,#34,#52;X3,X4](=O)([#6])[OX2H1]")
]


@contextlib.contextmanager
def acid_prefixes():
    """The acid-family prefix set (sulfonyl esters, boranyl groups ...) that substituents of a principal acid may use."""
    token = EXTENDED_PREFIXES.set(True)
    try:
        yield
    finally:
        EXTENDED_PREFIXES.reset(token)


_ELEMENT_RANK = {15: 0, 33: 1, 51: 2}
_ACID_CLASS_RANK = {(1, 2): 0, (0, 2): 1, (1, 1): 2, (0, 1): 3}
_NOT_AN_ACID = (1,)
_ACYL_PARTNERS = {7, 8, 9, 16, 17, 34, 35, 52, 53}


def acid_centre_rank(mol, atom):
    """Seniority key of a P, As or Sb atom (P-41, P-42.3): a phosphonic, phosphonous, phosphinic or phosphinous acid
    centre ranks by element, then acid type; any other centre (acid halide, amide, ester, chalcogenide, organyl
    derivative) is a class below every acid and shares one junior key."""
    if atom.GetAtomicNum() not in CENTER_STEMS or atom.GetFormalCharge() or atom.GetDegree() > 4:
        return _NOT_AN_ACID
    oxo = hydroxyl = organyl = 0
    for n in atom.GetNeighbors():
        double = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        if n.GetAtomicNum() == 8 and double and n.GetDegree() == 1:
            oxo += 1
        elif n.GetAtomicNum() == 8 and not double and n.GetDegree() == 1 and n.GetTotalNumHs() == 1:
            hydroxyl += 1
        elif n.GetAtomicNum() not in _ACYL_PARTNERS and not double:
            organyl += 1
        else:
            return _NOT_AN_ACID
    if hydroxyl + organyl + atom.GetTotalNumHs() != 3 or (oxo, hydroxyl) not in _ACID_CLASS_RANK:
        return _NOT_AN_ACID
    return 0, _ELEMENT_RANK[atom.GetAtomicNum()], _ACID_CLASS_RANK[(oxo, hydroxyl)]


def senior_acid_centres(mol, candidates):
    """The candidates whose acid is as senior as every other acid centre of the molecule (several when they tie)."""
    ranks = {a.GetIdx(): acid_centre_rank(mol, a) for a in mol.GetAtoms() if a.GetAtomicNum() in CENTER_STEMS}
    ranked = [c for c in candidates if ranks.get(c.GetIdx(), _NOT_AN_ACID) != _NOT_AN_ACID]
    best = min((ranks[c.GetIdx()] for c in ranked), default=None)
    if best is None or best != min(ranks.values()):
        return []
    return [c for c in ranked if ranks[c.GetIdx()] == best]


def require_phosphorus_acid_scope(mol, phosphorus):
    """A phosphorus acid centre that outranks every other group of the molecule (P-41): carboxylic and the
    sulfur-group acids are senior to phosphonic and phosphinic acids, which become phosphono prefixes."""
    in_nitro_group = set().union(*nitro_group_atoms(mol))
    for atom in mol.GetAtoms():
        if (atom.GetFormalCharge() != 0 and atom.GetIdx() not in in_nitro_group) or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if not senior_acid_centres(mol, [phosphorus]):
        raise UnsupportedStructure("the phosphorus centre is not the senior acid centre of the molecule")
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
    senior = senior_acid_centres(mol, phosphorus_atoms)
    if not senior:
        from ._phosphinic_acid import _phosphinic_acid_phosphorus_atoms, name_phosphinic_acid

        if senior_acid_centres(mol, _phosphinic_acid_phosphorus_atoms(mol)):
            return name_phosphinic_acid(mol)
    if len(senior) != 1 and len(phosphorus_atoms) != 1:
        raise UnsupportedStructure("more than one phosphonic acid group is not supported yet")
    phosphorus = senior[0] if len(senior) == 1 else phosphorus_atoms[0]

    require_phosphorus_acid_scope(mol, phosphorus)
    group_oxygens = {n.GetIdx() for n in phosphorus.GetNeighbors() if n.GetAtomicNum() == 8}

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    (root,) = [n for n in graph[phosphorus.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    with acid_prefixes(), cited_stereo_around(mol, phosphorus.GetIdx()):
        name, is_compound = name_branch(graph, root, phosphorus.GetIdx(), halogens, aromatic_atoms, mol=mol)
    prefix = enclose(name) if is_compound else name
    return f"{prefix}{centre_nuclides(phosphorus.GetIdx(), oxygens=True)}{CENTER_STEMS[phosphorus.GetAtomicNum()]}onic acid"
