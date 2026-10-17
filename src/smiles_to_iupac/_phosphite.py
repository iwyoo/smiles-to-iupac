"""Naming of trialkyl phosphite esters (P(OR)3, R groups identical or
mixed, P(III) -- no P=O, unlike `_phosphate.py`'s P(V) case), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-67.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  esters of mononuclear noncarbon oxoacids are named by citing the
  substituent (alkyl/aryl) group(s) as separate word(s), in alphanumeric
  order if more than one, followed by the acid's anion name. Confirmed
  worked example: `P(O-CH3)3` -> `trimethyl phosphite (PIN)`
  (the Blue Book) -- the exact same rule and
  citation style as the phosphate case, differing only in the anion word
  and the absence of a P=O bond.
- Word-assembly logic (grouping identical R names, multiplying-prefix +
  enclosure per group, alphanumeric word order) is shared verbatim with
  `_phosphate.py` via its exported `format_ester_words` -- see that
  module's own docstring for the full derivation.

Scope: a single phosphorus atom shaped like a phosphite ester -- no P=O
bond, three P-O-R single bonds, each R named via `name_branch` with the
groups that rank below the ester as prefixes (P-41), identical or mixed
freely. Explicitly out of scope (raise `UnsupportedStructure`): partial
("hydrogen") esters (a remaining P-OH, fewer than 3 ester oxygens), any
chalcogen-replacement analogue, and an acid, anhydride or other ester in R.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._phosphate import check_ester_organyl, cited_ester_groups, format_ester_entries

_PHOSPHORUS = 15


def _phosphite_phosphorus_atoms(mol):
    """Phosphorus atoms shaped like a fully-esterified phosphite: no
    double-bonded oxygen, exactly three single-bonded ester oxygens (each
    degree 2 -- P plus one R carbon; this excludes a hydroxyl oxygen,
    degree 1, i.e. a partial ester)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _PHOSPHORUS or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0 or atom.IsInRing():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(oxygens) != 3:
            continue
        if any(mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() != 1.0 for o in oxygens):
            continue
        if any(o.GetDegree() != 2 for o in oxygens):
            continue
        matches.append(atom)
    return matches


def has_phosphite_shape(mol) -> bool:
    return bool(_phosphite_phosphorus_atoms(mol))


def name_phosphite(mol) -> str:
    phosphorus_atoms = _phosphite_phosphorus_atoms(mol)
    if len(phosphorus_atoms) != 1:
        raise UnsupportedStructure("more than one phosphite group is not supported yet")
    (phosphorus,) = phosphorus_atoms

    group_oxygens = {n.GetIdx() for n in phosphorus.GetNeighbors() if n.GetAtomicNum() == 8}
    check_ester_organyl(mol, group_oxygens | {phosphorus.GetIdx()})

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    roots = []
    for oxygen_idx in group_oxygens:
        (root,) = [n for n in graph[oxygen_idx] if n != phosphorus.GetIdx()]
        roots.append((root, oxygen_idx))

    entries = cited_ester_groups(mol, graph, roots, halogens, aromatic_atoms)
    return format_ester_entries(entries) + " phosphite"
