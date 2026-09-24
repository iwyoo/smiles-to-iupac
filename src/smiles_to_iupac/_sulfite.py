"""Naming of dialkyl sulfite esters (S(=O)(OR)2, R groups identical or
mixed, one fewer double-bonded oxygen than `_sulfate.py`'s case), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-67.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  esters of mononuclear noncarbon oxoacids are named by citing the
  substituent (alkyl/aryl) group(s) as separate word(s), in alphanumeric
  order if more than one, followed by the acid's anion name -- mirrors
  the phosphate/phosphite relationship, same central-atom/citation-style
  pair as `_sulfate.py`, differing only in the anion word ("sulfite")
  and the absence of a second S=O bond.
- Word-assembly logic is shared verbatim with `_phosphate.py` via its
  exported `format_ester_words`.

Scope: a single sulfur atom shaped like a sulfite ester -- one S=O
double bond, two S-O-R single bonds, each R named via `name_branch` (a
plain alkyl chain, a branched chain, or a plain benzene ring, and their
halogenated variants, exactly like `_sulfate.py`'s R), identical or
mixed freely. Explicitly out of scope (raise `UnsupportedStructure`):
partial ("hydrogen") esters (a remaining S-OH), any chalcogen-
replacement analogue, and any other heteroatom.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._phosphate import format_ester_words
from ._substituents import name_branch

_SULFUR = 16


def _sulfite_sulfur_atoms(mol):
    """Sulfur atoms shaped like a fully-esterified sulfite: bonded to
    exactly one double-bonded (terminal) oxygen and exactly two single-
    bonded ester oxygens (each degree 2 -- S plus one R carbon; this
    excludes a hydroxyl oxygen, degree 1, i.e. a partial ester)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _SULFUR or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(oxygens) != 3:
            continue
        double_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(ester_os) != 2:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        if any(o.GetDegree() != 2 for o in ester_os):
            continue
        matches.append(atom)
    return matches


def has_sulfite_shape(mol) -> bool:
    return bool(_sulfite_sulfur_atoms(mol))


def name_sulfite(mol) -> str:
    sulfur_atoms = _sulfite_sulfur_atoms(mol)
    if len(sulfur_atoms) != 1:
        raise UnsupportedStructure("more than one sulfite group is not supported yet")
    (sulfur,) = sulfur_atoms

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _SULFUR and atom.GetIdx() != sulfur.GetIdx():
            raise UnsupportedStructure("more than one sulfur atom is not supported yet")
        if atomic_num not in (1, 6, 8, _SULFUR, *HALOGEN_PREFIXES):
            raise UnsupportedStructure(
                "heteroatoms other than the sulfite's own sulfur/"
                "oxygens and a halogen substituent are not supported yet"
            )

    group_oxygens = {n.GetIdx() for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetIdx() not in group_oxygens:
            raise UnsupportedStructure(
                "an oxygen atom not part of the sulfite's own "
                "S(=O)(OR)2 group is out of scope for this module"
            )

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    double_os = {
        n.GetIdx()
        for n in sulfur.GetNeighbors()
        if n.GetAtomicNum() == 8
        and mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    }
    ester_oxygens = [idx for idx in group_oxygens if idx not in double_os]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    roots = []
    for oxygen_idx in ester_oxygens:
        (root,) = [n for n in graph[oxygen_idx] if n != sulfur.GetIdx()]
        roots.append((root, oxygen_idx))

    names = [name_branch(graph, root, coming_from, halogens, aromatic_atoms, mol=mol)[0] for root, coming_from in roots]
    return format_ester_words(names) + " sulfite"
