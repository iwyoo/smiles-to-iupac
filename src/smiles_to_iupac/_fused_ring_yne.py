"""Ring triple bonds in a fused or fused-bridged ring system (P-31.2.1, P-31.2.4.1): the mancude or hydro parent with the
prefix 'didehydro', 'tetradehydro', ... as in '2,3-didehydro-1H-indene' or '3,4-didehydro-1,2-dihydroquinoline'. The
locants are read from two probes that turn each triple bond into a double bond: one carrying iodo prefixes (hydro and
dehydro prefixes together get the lowest locants) and one carrying hydroxy suffixes (a tie goes to the lowest dehydro
locants, P-31.1.4.3.4)."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure

_HYDRO_GROUP = r"(?:(?:\d+[a-z]?,)*\d+[a-z]?-(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?hydro-?)"
_INDICATED_GROUP = r"(?:\d+[a-z]?H-)"
_HYDRO = re.compile(rf"^({_HYDRO_GROUP})")
_REFUSED = ("spiro", "bicyclo", "tricyclo", "phane", "iodo", "hydroxy", "diol")
_BRIDGE = re.compile(r"\d+,\d+-[a-z]+ano")


_MULTIPLIER = {1: "di", 2: "tetra", 3: "hexa"}


def _triple_bonds(mol):
    bonds = [b for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 3.0]
    if not bonds or len(bonds) > len(_MULTIPLIER) or not all(b.IsInRing() for b in bonds):
        return None
    if any(a.GetAtomicNum() != 6 for b in bonds for a in (b.GetBeginAtom(), b.GetEndAtom())):
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    return bonds if len(Chem.GetMolFrags(mol)) == 1 else None


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
    bonds = _triple_bonds(mol)
    if bonds is None:
        return False
    if all(_ortho_fused_system(mol, b) for b in bonds):
        return True
    word = _parent_word(mol, bonds)
    return word is not None and not any(refused in word for refused in _REFUSED) and bool(_BRIDGE.search(word))


def _name_of(editable):
    from .core import smiles_to_iupac

    probe = editable.GetMol()
    try:
        Chem.SanitizeMol(probe)
        return smiles_to_iupac(Chem.MolToSmiles(probe))
    except (UnsupportedStructure, ValueError):
        return None


def _probe_name(mol, bonds, atomic_number):
    editable = Chem.RWMol(mol)
    for bond in bonds:
        editable.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()).SetBondType(Chem.BondType.DOUBLE)
        for atom in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()):
            added = editable.AddAtom(Chem.Atom(atomic_number))
            editable.AddBond(atom, added, Chem.BondType.SINGLE)
    return _name_of(editable)


def _parent_word(mol, bonds):
    """The parent hydride without hydro prefixes and indicated hydrogen, read from the ring atoms alone."""
    editable = Chem.RWMol(mol)
    for bond in bonds:
        editable.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()).SetBondType(Chem.BondType.DOUBLE)
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    for index in sorted(set(range(mol.GetNumAtoms())) - ring_atoms, reverse=True):
        editable.RemoveAtom(index)
    name = _name_of(editable)
    if name is None:
        return None
    return re.sub(rf"^{_HYDRO_GROUP}?{_INDICATED_GROUP}*", "", name)


def _split(rest, parent_word):
    """(prefixes, hydro/indicated hydrogen with the parent and any suffix, length of the hydro and indicated part) of a
    probe name without its probe marker."""
    stem = parent_word[:-1] if parent_word.endswith("e") else parent_word
    match = re.search(rf"(?P<hydro>{_HYDRO_GROUP}?{_INDICATED_GROUP}*)(?:{re.escape(parent_word)}|{re.escape(stem)}(?=-))", rest)
    if match is None:
        return None
    return rest[: match.start()], rest[match.start():], len(match.group("hydro"))


def _candidate(dehydro, rest, parent_word):
    """(sort key, name) of the didehydro name built from a probe name, None when the probe name has another form."""
    if any(word in rest for word in _REFUSED):
        return None
    parts = _split(rest, parent_word)
    if parts is None:
        return None
    head, tail, block = parts
    hydro_match = _HYDRO.match(tail)
    hydro = [int(n) for n in re.findall(r"\d+", hydro_match.group(1).split("hydro")[0])] if hydro_match else []
    indicated = [int(n) for n in re.findall(r"(\d+)[a-z]?H-", tail[:block])]
    prefix = ",".join(map(str, sorted(dehydro))) + "-" + _MULTIPLIER[len(dehydro) // 2] + "dehydro"
    text = head + ("-" if head and not head.endswith("-") else "") + prefix + ("-" if tail[0].isdigit() else "") + tail
    return (indicated, sorted(hydro + dehydro), dehydro), text


def name_fused_ring_yne(mol) -> str:
    bonds = _triple_bonds(mol)
    if bonds is None:
        raise UnsupportedStructure("not a ring triple bond of a fused carbocycle or heterocycle")
    parent_word = _parent_word(mol, bonds)
    if parent_word is None:
        raise UnsupportedStructure("the ring system of this ring triple bond has no plain parent hydride name")
    multiplier = _MULTIPLIER[len(bonds)]
    found = []
    iodo = _probe_name(mol, bonds, 53)
    match = re.search(rf"(?:^|-)((?:\d+,)+\d+)-{multiplier}iodo-?", iodo) if iodo else None
    if match:
        before, after = iodo[: match.start()], iodo[match.end():]
        rest = before + ("-" if before and after[:1].isdigit() else "") + after
        candidate = _candidate([int(n) for n in match.group(1).split(",")], rest, parent_word)
        if candidate:
            found.append(candidate)
    # a hydroxy probe would outrank an amino or other suffix group of the molecule and change the parent
    diol = _probe_name(mol, bonds, 8) if all(a.GetAtomicNum() in (6, 9, 17, 35) or a.IsInRing() for a in mol.GetAtoms()) else None
    match = re.search(rf"-((?:\d+,)+\d+)-{multiplier}ol$", diol) if diol else None
    if match:
        candidate = _candidate([int(n) for n in match.group(1).split(",")], diol[: match.start()], parent_word)
        if candidate:
            found.append(candidate)
    if not found:
        raise UnsupportedStructure("the dehydro prefixes of this fused ring system could not be placed")
    return min(found, key=lambda c: c[0])[1]
