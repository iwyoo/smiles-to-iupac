"""Naming of nitrite esters (R-O-N=O, a single R group), per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-67.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  esters of mononuclear noncarbon oxoacids are named by citing the
  substituent (alkyl/aryl) group as a separate word, followed by the
  acid's anion name -- mirrors `_nitrate_ester.py`'s own citation style
  and module shape (a separate module, not a shared detector
  parameterized by double-bond-O count, mirroring the established
  sulfate/sulfite precedent), differing only in the anion word
  ("nitrite") and one fewer oxygen.
- Word-assembly logic is shared verbatim with `_phosphate.py` via its
  exported `format_ester_words` -- exercised here only in its simplest
  (single word, no multiplying prefix) form, same as `_nitrate_ester.py`.
- Unlike `_nitrate_ester.py`'s nitrogen, there is no charged/neutral
  dual-representation ambiguity here -- nitrite's N=O has no analogous
  resonance form needing separate handling, so the nitrogen's formal
  charge is always 0.

Scope: a single nitrogen atom shaped like a nitrite ester -- one N-O-R
single bond, one N=O double bond, R named via `name_branch` (a plain
alkyl chain, a branched chain, or a plain benzene ring, and their
halogenated variants, exactly like `_nitrate_ester.py`'s R). Explicitly
out of scope (raise `UnsupportedStructure`): any chalcogen-replacement
analogue, any other heteroatom, more than one nitrite group.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._phosphate import format_ester_words
from ._substituents import name_branch

_NITROGEN = 7
_OXYGEN = 8


def _nitrite_ester_nitrogen_atoms(mol):
    """Nitrogen atoms shaped like a nitrite ester: bonded to zero
    carbons, one single-bonded ester oxygen (degree 2 -- N plus one R
    carbon) and one double-bonded terminal oxygen (degree 1), formal
    charge 0."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _NITROGEN or atom.GetDegree() != 2 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        if any(n.GetAtomicNum() == 6 for n in neighbors):
            continue
        oxygens = [n for n in neighbors if n.GetAtomicNum() == _OXYGEN]
        if len(oxygens) != 2:
            continue
        ester_os = [
            o
            for o in oxygens
            if o.GetDegree() == 2 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        terminal_os = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(ester_os) != 1 or len(terminal_os) != 1:
            continue
        matches.append(atom)
    return matches


def has_nitrite_ester_shape(mol) -> bool:
    return bool(_nitrite_ester_nitrogen_atoms(mol))


def name_nitrite_ester(mol) -> str:
    nitrogen_atoms = _nitrite_ester_nitrogen_atoms(mol)
    if len(nitrogen_atoms) != 1:
        raise UnsupportedStructure("more than one nitrite group is not supported yet")
    (nitrogen,) = nitrogen_atoms

    group_atom_idxs = {nitrogen.GetIdx()} | {n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == _OXYGEN}
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:
            raise UnsupportedStructure("isotopically modified atoms are not supported yet")
        if atom.GetIdx() not in group_atom_idxs and atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("charged atoms are not supported yet")
        if atom.GetIdx() == nitrogen.GetIdx():
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _NITROGEN:
            raise UnsupportedStructure("more than one nitrogen atom is not supported yet")
        if atomic_num not in (1, 6, _OXYGEN, *HALOGEN_PREFIXES):
            raise UnsupportedStructure(
                "heteroatoms other than the nitrite's own nitrogen/"
                "oxygens and a halogen substituent are not supported yet"
            )

    group_oxygens = {n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == _OXYGEN}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == _OXYGEN and atom.GetIdx() not in group_oxygens:
            raise UnsupportedStructure(
                "an oxygen atom not part of the nitrite's own "
                "N(=O)OR group is out of scope for this module"
            )

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (ester_oxygen_idx,) = [
        idx
        for idx in group_oxygens
        if mol.GetAtomWithIdx(idx).GetDegree() == 2
        and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), idx).GetBondTypeAsDouble() == 1.0
    ]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    (root,) = [n for n in graph[ester_oxygen_idx] if n != nitrogen.GetIdx()]
    name, _ = name_branch(graph, root, ester_oxygen_idx, halogens, aromatic_atoms, mol=mol)
    return format_ester_words([name]) + " nitrite"
