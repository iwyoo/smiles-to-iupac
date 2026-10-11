"""A ring triple bond in an ortho-fused carbocycle (P-31.2.1, P-31.2.4.1): the mancude or hydro parent with the prefix
'didehydro', '2,3-didehydro-1H-indene', '1,2-didehydro-3,4-dihydronaphthalene'. The locants are read from two probes
that turn the triple bond into a double bond: one carrying iodo prefixes (hydro and dehydro prefixes together get the
lowest locants) and one carrying hydroxy suffixes (a tie goes to the lowest dehydro locants, P-31.1.4.3.4)."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure

_HYDRO_GROUP = r"(?:(?:\d+[a-z]?,)*\d+[a-z]?-(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?hydro-?)"
_INDICATED_GROUP = r"(?:\d+[a-z]?H-)"
_HYDRO = re.compile(rf"^({_HYDRO_GROUP})")
_REFUSED = ("spiro", "bicyclo", "tricyclo", "phane", "iodo", "hydroxy", "diol")


def _triple_bond(mol):
    bonds = [b for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 3.0]
    if len(bonds) != 1 or not bonds[0].IsInRing():
        return None
    if any(a.GetAtomicNum() != 6 or a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    return bonds[0] if len(Chem.GetMolFrags(mol)) == 1 else None


def _ortho_fused_system(mol, bond):
    """Whether the rings joined to the ring of `bond` share bonds only (no spiro atom, no bridge)."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    start = [r for r in rings if bond.GetBeginAtomIdx() in r and bond.GetEndAtomIdx() in r]
    if not start:
        return False
    system, grown = [start[0]], True
    while grown:
        grown = False
        for r in rings:
            if r not in system and any(r & s for s in system):
                system.append(r)
                grown = True
    if len(system) < 2:
        return False
    return all(len(a & b) in (0, 2) for i, a in enumerate(system) for b in system[:i])


def has_fused_ring_yne_shape(mol) -> bool:
    bond = _triple_bond(mol)
    return bond is not None and _ortho_fused_system(mol, bond)


def _name_of(editable):
    from .core import smiles_to_iupac

    probe = editable.GetMol()
    try:
        Chem.SanitizeMol(probe)
        return smiles_to_iupac(Chem.MolToSmiles(probe))
    except (UnsupportedStructure, ValueError):
        return None


def _probe_name(mol, bond, atomic_number):
    editable = Chem.RWMol(mol)
    editable.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()).SetBondType(Chem.BondType.DOUBLE)
    for atom in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()):
        added = editable.AddAtom(Chem.Atom(atomic_number))
        editable.AddBond(atom, added, Chem.BondType.SINGLE)
    return _name_of(editable)


def _parent_word(mol, bond):
    """The parent hydride without hydro prefixes and indicated hydrogen, read from the ring atoms alone."""
    editable = Chem.RWMol(mol)
    editable.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()).SetBondType(Chem.BondType.DOUBLE)
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    if any(
        b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds() if (b.GetBeginAtomIdx() in ring_atoms) != (b.GetEndAtomIdx() in ring_atoms)
    ):
        return None
    for index in sorted(set(range(mol.GetNumAtoms())) - ring_atoms, reverse=True):
        editable.RemoveAtom(index)
    name = _name_of(editable)
    if name is None:
        return None
    return re.sub(rf"^{_HYDRO_GROUP}?{_INDICATED_GROUP}*", "", name)


def _split(rest, parent_word):
    """(prefixes, hydro/indicated hydrogen and parent) of a probe name without its probe marker."""
    match = re.search(rf"{_HYDRO_GROUP}?{_INDICATED_GROUP}*{re.escape(parent_word)}$", rest)
    return (rest[: match.start()], rest[match.start():]) if match else None


def _candidate(dehydro, rest, parent_word):
    """(sort key, name) of the didehydro name built from a probe name, None when the probe name has another form."""
    if any(word in rest for word in _REFUSED):
        return None
    parts = _split(rest, parent_word)
    if parts is None:
        return None
    head, tail = parts
    hydro_match = _HYDRO.match(tail)
    hydro = [int(n) for n in re.findall(r"\d+", hydro_match.group(1).split("hydro")[0])] if hydro_match else []
    indicated = [int(n) for n in re.findall(rf"(\d+)[a-z]?H-", tail[: len(tail) - len(parent_word)])]
    prefix = f"{dehydro[0]},{dehydro[1]}-didehydro"
    text = head + ("-" if head and not head.endswith("-") else "") + prefix + ("-" if tail[0].isdigit() else "") + tail
    return (indicated, sorted(hydro + dehydro), dehydro), text


def name_fused_ring_yne(mol) -> str:
    bond = _triple_bond(mol)
    if bond is None or not _ortho_fused_system(mol, bond):
        raise UnsupportedStructure("not a ring triple bond of an ortho-fused carbocycle")
    parent_word = _parent_word(mol, bond)
    if parent_word is None:
        raise UnsupportedStructure("the ring system of this ring triple bond has no plain parent hydride name")
    found = []
    iodo = _probe_name(mol, bond, 53)
    match = re.search(r"(?:^|-)(\d+),(\d+)-diiodo-?", iodo) if iodo else None
    if match:
        before, after = iodo[: match.start()], iodo[match.end():]
        rest = before + ("-" if before and after[:1].isdigit() else "") + after
        candidate = _candidate([int(match.group(1)), int(match.group(2))], rest, parent_word)
        if candidate:
            found.append(candidate)
    diol = _probe_name(mol, bond, 8)
    match = re.search(r"-(\d+),(\d+)-diol$", diol) if diol else None
    if match:
        candidate = _candidate([int(match.group(1)), int(match.group(2))], diol[: match.start()], parent_word)
        if candidate:
            found.append(candidate)
    if not found:
        raise UnsupportedStructure("the didehydro prefixes of this fused ring system could not be placed")
    return min(found, key=lambda c: c[0])[1]
