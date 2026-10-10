"""Two identical nitrogen-bearing parents joined through their nitrogen atoms by an alkanediyl group (P-62.2.5.2):
N1,N1'-methylenedi(ethane-1,2-diamine), N',N'''-(ethane-1,2-diyl)diacetohydrazide.
Amine parents are only joined multiplicatively through methylene; a longer chain between two amine nitrogens
is itself a diamine parent (P-62.2.4.1.2), and so is a methylene between two monoamine nitrogens (methanediamine)."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._multiplicative_text import unit_phrase
from ._numerals import alkane_name

_MARKER = "tert-butyl"
_MARKER_LOCANT = re.compile(r"((?:N['′]*\d*,)*N['′]*\d*)-(?:\d+[a-z]*-)?" + _MARKER)


def _is_imine_nitrogen(atom):
    return (
        atom.GetAtomicNum() == 7
        and not atom.IsInRing()
        and atom.GetDegree() == 2
        and not atom.GetTotalNumHs()
        and sorted(b.GetBondTypeAsDouble() for b in atom.GetBonds()) == [1.0, 2.0]
        and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds())
    )


def _is_attach_nitrogen(atom):
    return _is_imine_nitrogen(atom) or (
        atom.GetAtomicNum() == 7 and not atom.IsInRing() and atom.GetTotalNumHs() == 1 and atom.GetDegree() == 2
    )


def _is_methylene(atom):
    return atom.GetAtomicNum() == 6 and not atom.IsInRing() and atom.GetTotalNumHs() == 2 and atom.GetDegree() == 2


def _imine_carbon_with_ring_substituent(smiles):
    unit = Chem.MolFromSmiles(smiles)
    return any(
        a.GetAtomicNum() == 6 and not a.IsInRing() and any(n.IsInRing() for n in a.GetNeighbors())
        for a in unit.GetAtoms()
        if any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(a).GetAtomicNum() == 7 for b in a.GetBonds())
    )


def _arms(mol):
    for start in mol.GetAtoms():
        if not _is_attach_nitrogen(start):
            continue
        for first in start.GetNeighbors():
            if mol.GetBondBetweenAtoms(start.GetIdx(), first.GetIdx()).GetBondTypeAsDouble() != 1.0:
                continue
            path, previous, current = [], start, first
            while _is_methylene(current):
                path.append(current)
                following = next(n for n in current.GetNeighbors() if n.GetIdx() != previous.GetIdx())
                previous, current = current, following
            if path and _is_attach_nitrogen(current) and start.GetIdx() < current.GetIdx():
                yield path, [start, current]


def _side(mol, blocked, nitrogen):
    seen, stack = {nitrogen.GetIdx()}, [nitrogen.GetIdx()]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetIdx() not in seen and n.GetIdx() not in blocked:
                seen.add(n.GetIdx())
                stack.append(n.GetIdx())
    return seen


def _fragment(mol, atoms, attach, marked):
    editable = Chem.RWMol(mol)
    extra = []
    if marked:
        tail = editable.AddAtom(Chem.Atom(6))
        editable.AddBond(attach, tail, Chem.BondType.SINGLE)
        extra.append(tail)
        for _ in range(3):
            methyl = editable.AddAtom(Chem.Atom(6))
            editable.AddBond(tail, methyl, Chem.BondType.SINGLE)
            extra.append(methyl)
    keep = set(atoms) | set(extra)
    for idx in sorted(set(range(editable.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(idx)
    fragment = editable.GetMol()
    Chem.SanitizeMol(fragment)
    return Chem.MolToSmiles(fragment)


def _nth_locant(token, index, nitrogens_per_unit):
    letters = token[1:]
    digits = re.sub(r"['′]", "", letters)
    if digits:
        return f"N{digits}{chr(39) * index}"
    return "N" + chr(39) * (len(letters) + index * nitrogens_per_unit)


def _has_stereo(mol):
    return any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )


def name_nitrogen_methylene_multiplicative(mol):
    from .core import smiles_to_iupac

    if len(Chem.GetMolFrags(mol)) != 1 or _has_stereo(mol):
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    for path, pair in _arms(mol):
        blocked = {a.GetIdx() for a in path}
        sides = [_side(mol, blocked, n) for n in pair]
        if sides[0] & sides[1] or len(sides[0]) + len(sides[1]) + len(path) != mol.GetNumAtoms():
            continue
        if any(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for s in sides for i in s):
            continue
        plain = [_fragment(mol, s, n.GetIdx(), False) for s, n in zip(sides, pair)]
        if plain[0] != plain[1]:
            continue
        nitrogens = sum(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in sides[0])
        try:
            unit_name = smiles_to_iupac(plain[0])
            marked_name = (
                "" if nitrogens == 1 else smiles_to_iupac(_fragment(mol, sides[0], pair[0].GetIdx(), True))
            )
        except UnsupportedStructure:
            continue
        found = _MARKER_LOCANT.search(marked_name)
        if _MARKER in unit_name or (nitrogens > 1 and found is None):
            continue
        if len(path) > 1 and not unit_name.endswith(("hydrazide", "imine")):
            continue
        if nitrogens == 1 and len(path) == 1 and unit_name.endswith(("amine", "aniline")):
            continue
        token = "N" if nitrogens == 1 else found.group(1).split(",")[-1]
        locants = ",".join(_nth_locant(token, k, nitrogens) for k in range(2))
        central = "methylene" if len(path) == 1 else f"({alkane_name(len(path))}-1,{len(path)}-diyl)"
        if unit_name.endswith("imine") and _imine_carbon_with_ring_substituent(plain[0]):
            return f"{locants}-{central}bis({unit_name})"
        return f"{locants}-{central}di{unit_phrase(unit_name, '')}"
    return None
