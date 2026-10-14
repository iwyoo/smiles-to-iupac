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
single bond, one N=O double bond, R named via `name_branch` with the groups
that rank below the ester as prefixes, as for `_nitrate_ester.py`'s R.
Out of scope (raise `UnsupportedStructure`): any chalcogen-replacement
analogue, R holding an acid or another ester. Several nitrite groups on one organyl group are named like the
acetates of the same polyol (P-65.6.3.2.2).
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._phosphate import check_ester_organyl, format_ester_words, name_polyester_through_acetates
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
    if len(nitrogen_atoms) > 1:
        group_atoms = {
            n.GetIdx() for nitrogen in nitrogen_atoms for n in [nitrogen, *nitrogen.GetNeighbors()] if n.GetAtomicNum() in (7, 8)
        }
        check_ester_organyl(mol, group_atoms)
        return name_polyester_through_acetates(mol, [n.GetIdx() for n in nitrogen_atoms], "nitrite", "nitrosooxy")
    (nitrogen,) = nitrogen_atoms

    group_atom_idxs = {nitrogen.GetIdx()} | {n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == _OXYGEN}
    check_ester_organyl(mol, group_atom_idxs)
    group_oxygens = {n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == _OXYGEN}

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
