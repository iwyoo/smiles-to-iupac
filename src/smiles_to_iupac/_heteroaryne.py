"""Didehydro derivatives of mancude heterocycles, the heteroarynes (P-31.1.4.2.4, P-44.4.1.10.2): a ring C#C bond that
becomes aromatic when read as C=C is named from the parent ring with the prefix 'didehydro', '3,4-didehydropyridine'.
The locants are read from the diiodo derivative of the parent."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure


def _triple_bond(mol):
    bonds = [b for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 3.0]
    if len(bonds) != 1 or not bonds[0].IsInRing():
        return None
    bond = bonds[0]
    if any(a.GetAtomicNum() != 6 for a in (bond.GetBeginAtom(), bond.GetEndAtom())):
        return None
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    return bond


def has_heteroaryne_shape(mol) -> bool:
    bond = _triple_bond(mol)
    if bond is None or all(a.GetAtomicNum() == 6 for a in mol.GetAtoms()):
        return False
    return _parent_probe(mol, bond) is not None


def _parent_probe(mol, bond):
    editable = Chem.RWMol(mol)
    target = editable.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    target.SetBondType(Chem.BondType.DOUBLE)
    for atom in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()):
        iodine = editable.AddAtom(Chem.Atom(53))
        editable.AddBond(atom, iodine, Chem.BondType.SINGLE)
    probe = editable.GetMol()
    try:
        Chem.SanitizeMol(probe)
    except Exception:
        return None
    ring_atoms = {a for ring in probe.GetRingInfo().AtomRings() for a in ring}
    if not all(probe.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    return probe


def name_heteroaryne(mol) -> str:
    from .core import smiles_to_iupac

    bond = _triple_bond(mol)
    probe = _parent_probe(mol, bond) if bond is not None else None
    if probe is None:
        raise UnsupportedStructure("not a didehydro derivative of an aromatic ring")
    name = smiles_to_iupac(Chem.MolToSmiles(probe))
    match = re.fullmatch(r"(\d+[a-z]?),(\d+[a-z]?)-diiodo(.+)", name)
    if match is None:
        raise UnsupportedStructure("the locants of the didehydro prefix could not be read")
    return f"{match.group(1)},{match.group(2)}-didehydro{match.group(3)}"
