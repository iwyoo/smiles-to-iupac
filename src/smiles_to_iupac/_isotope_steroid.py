"""Isotopically modified steroid parent hydrides (P-82.4.1, P-82.4.2, P-82.5.2): the unmodified parent is named, the
descriptor of the labelled positions is cited with their steroid locants in front of the parent hydride name, and a
centre that the nuclide makes stereogenic is cited first: (24R)-5α-(24-2H1)cholestane, 5α-(17-2H)pregnane."""

import re

from rdkit import Chem

from ._axial_stereo import _merge
from ._common import UnsupportedStructure, assign_cip_labels
from ._isotope_labels import descriptor, modification_key, split_isotopes
from ._steroid_named import _core_atoms, _embeddings

_PARENT = re.compile(r"(gon|estr|androst|pregn|cholest|chol|ergost|campest|poriferast|stigmast)ane$")
_TEMPLATE = {"campest": "ergost", "poriferast": "stigmast"}


def has_isotope_steroid_shape(mol) -> bool:
    return (
        mol.GetRingInfo().NumRings() == 4
        and any(a.GetIsotope() for a in mol.GetAtoms())
        and all(a.GetAtomicNum() in (1, 6) for a in mol.GetAtoms())
    )


def _positions(clean, stem, labels):
    graph = {a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors()} for a in clean.GetAtoms()}
    ring_atoms = _core_atoms(clean, graph)
    if ring_atoms is None:
        raise UnsupportedStructure("this isotopically modified ring system is not a steroid skeleton")
    best = None
    for mapping in _embeddings(_TEMPLATE.get(stem, stem), graph, ring_atoms):
        position_of = {atom: position for position, atom in mapping.items()}
        if any(atom not in position_of for atom in labels):
            continue
        key = modification_key(labels, position_of)
        if best is None or key < best[0]:
            best = (key, position_of)
    if best is None:
        raise UnsupportedStructure("an isotopically modified atom outside the steroid parent hydride")
    return best[1]


def _new_centres(mol, clean, index, position_of):
    """(locant, R/S) of the centres that are stereogenic only through the nuclide."""
    assign_cip_labels(mol)
    specified = Chem.StereoSpecified.Specified
    tetrahedral = Chem.StereoType.Atom_Tetrahedral
    inherent = {e.centeredOn for e in Chem.FindPotentialStereo(clean) if e.type == tetrahedral and e.specified == specified}
    found = []
    for element in Chem.FindPotentialStereo(mol):
        if element.type != tetrahedral or element.specified != specified or index[element.centeredOn] in inherent:
            continue
        atom = mol.GetAtomWithIdx(element.centeredOn)
        code = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        if code not in ("R", "S") or index[element.centeredOn] not in position_of:
            raise UnsupportedStructure("the stereocentre made by this nuclide is not cited by a supported name")
        found.append((position_of[index[element.centeredOn]], code))
    return found


def name_isotope_steroid(mol) -> str:
    from .core import smiles_to_iupac

    split = split_isotopes(mol)
    if split is None:
        raise UnsupportedStructure("no isotopic label")
    clean, labels, index = split
    name = smiles_to_iupac(Chem.MolToSmiles(clean))
    parent = _PARENT.search(name)
    if parent is None:
        raise UnsupportedStructure("the unmodified steroid is not named as a bare parent hydride")
    position_of = _positions(clean, parent.group(1), labels)
    capacity = {atom: clean.GetAtomWithIdx(atom).GetTotalNumHs() for atom in labels}
    text = descriptor(labels, position_of, False, capacity=capacity)
    name = name[: parent.start()] + text + name[parent.start() :]
    for locant, code in sorted(_new_centres(mol, clean, index, position_of)):
        name = _merge(name, f"{locant}{code}")
    return name
