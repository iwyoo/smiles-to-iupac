"""Polycations with identical ylium centres on characteristic groups of one skeleton (P-73.5.1.2): carboxamidylium
centres (-C(=O)-NH+), sulfonylium centres (-SO2+, cited as dioxo-lambda6-sulfanylium) and disulfanylium centres on acyl
groups (-C(=O)-S-S+). The skeleton is named from a surrogate with the same attachments: the neutral polyamide, the
diacetate of the diol, or the diacyl dichloride."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._numerals import multiplying_prefix

_ALLOWED = (1, 6, 7, 8, 9, 17, 35, 53)


def _carbonyl_carbon(atom):
    return next(
        (
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 6
            and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(n).GetAtomicNum() == 8 for b in n.GetBonds())
        ),
        None,
    )


def _groups(mol):
    """(kind, [(anchor, members)]) of the identical ylium groups of `mol`, else None."""
    centers = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(centers) < 2 or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() != 1 or a.IsInRing() or a.GetIsotope() for a in centers):
        return None
    found = []
    for atom in centers:
        z = atom.GetAtomicNum()
        neighbors = list(atom.GetNeighbors())
        bonds = {n.GetIdx(): mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in neighbors}
        if z == 7 and len(neighbors) == 1 and atom.GetTotalNumHs() == 1 and bonds[neighbors[0].GetIdx()] == 1.0:
            anchor = _carbonyl_carbon(atom)
            if anchor is None:
                return None
            found.append(("amidylium", anchor.GetIdx(), {atom.GetIdx()}))
        elif z == 16 and len(neighbors) == 3 and not atom.GetTotalNumHs():
            oxo = [n for n in neighbors if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and bonds[n.GetIdx()] == 2.0]
            rest = [n for n in neighbors if n not in oxo]
            if len(oxo) != 2 or len(rest) != 1 or rest[0].GetAtomicNum() != 6 or bonds[rest[0].GetIdx()] != 1.0:
                return None
            found.append(("sulfonylium", rest[0].GetIdx(), {atom.GetIdx(), *(o.GetIdx() for o in oxo)}))
        elif z == 16 and len(neighbors) == 1 and not atom.GetTotalNumHs() and neighbors[0].GetAtomicNum() == 16:
            host = neighbors[0]
            carbon = [n for n in host.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
            if host.GetDegree() != 2 or len(carbon) != 1 or carbon[0].GetAtomicNum() != 6:
                return None
            if not any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(carbon[0]).GetAtomicNum() == 8 for b in carbon[0].GetBonds()):
                return None
            found.append(("disulfanylium", carbon[0].GetIdx(), {atom.GetIdx(), host.GetIdx()}))
        else:
            return None
    if len({kind for kind, _, _ in found}) != 1:
        return None
    members = set().union(*(m for _, _, m in found))
    if any(a.GetAtomicNum() not in _ALLOWED and a.GetIdx() not in members for a in mol.GetAtoms()):
        return None
    if any(a.GetNumRadicalElectrons() and a.GetIdx() not in {c.GetIdx() for c in centers} for a in mol.GetAtoms()):
        return None
    return found[0][0], [(anchor, group) for _, anchor, group in found]


def has_group_polycation_shape(mol) -> bool:
    return _groups(mol) is not None


def _surrogate(mol, kind, groups):
    editable = Chem.RWMol(mol)
    removed = set()
    for anchor, members in groups:
        removed |= members
        if kind == "amidylium":
            (nitrogen,) = members
            atom = editable.GetAtomWithIdx(nitrogen)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(2)
            atom.SetNumRadicalElectrons(0)
        elif kind == "sulfonylium":
            sulfur = next(i for i in members if mol.GetAtomWithIdx(i).GetAtomicNum() == 16)
            atom = editable.GetAtomWithIdx(sulfur)
            atom.SetAtomicNum(8)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
            atom.SetNumRadicalElectrons(0)
            acyl = editable.AddAtom(Chem.Atom(6))
            oxo = editable.AddAtom(Chem.Atom(8))
            methyl = editable.AddAtom(Chem.Atom(6))
            editable.AddBond(sulfur, acyl, Chem.BondType.SINGLE)
            editable.AddBond(acyl, oxo, Chem.BondType.DOUBLE)
            editable.AddBond(acyl, methyl, Chem.BondType.SINGLE)
        else:
            end = next(i for i in members if mol.GetAtomWithIdx(i).GetFormalCharge())
            atom = editable.GetAtomWithIdx(end)
            atom.SetAtomicNum(17)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
            atom.SetNumRadicalElectrons(0)
            host = next(i for i in members if i != end)
            editable.RemoveBond(anchor, host)
            editable.AddBond(anchor, end, Chem.BondType.SINGLE)
            removed.add(host)
    if kind == "sulfonylium":
        removed = {i for i in removed if mol.GetAtomWithIdx(i).GetAtomicNum() == 8 and mol.GetAtomWithIdx(i).GetDegree() == 1}
    elif kind == "disulfanylium":
        removed = {i for i in removed if mol.GetAtomWithIdx(i).GetAtomicNum() == 16 and not mol.GetAtomWithIdx(i).GetFormalCharge()}
    else:
        removed = set()
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    surrogate = editable.GetMol()
    Chem.SanitizeMol(surrogate)
    return surrogate


def name_group_polycation(mol) -> str:
    from .core import smiles_to_iupac

    shape = _groups(mol)
    if shape is None:
        raise UnsupportedStructure("not a polycation of identical ylium groups")
    kind, groups = shape
    count = len(groups)
    name = smiles_to_iupac(Chem.MolToSmiles(_surrogate(mol, kind, groups)))
    if kind == "amidylium":
        match = re.search(r"(di|tri|tetra)(carbox)?amide$", name)
        if match is None:
            raise UnsupportedStructure("the skeleton of this polyamidylium has no polyamide name")
        multiplier = {"di": "bis", "tri": "tris", "tetra": "tetrakis"}[match.group(1)]
        return f"{name[: match.start()]}{multiplier}({match.group(2) or ''}amidylium)"
    if kind == "sulfonylium":
        tail = " " + multiplying_prefix(count) + "acetate"
        word = "(dioxo-λ6-sulfanylium)"
    else:
        tail = " " + multiplying_prefix(count) + "chloride"
        word = "(disulfanylium)"
    if not name.endswith(tail):
        raise UnsupportedStructure("the skeleton of this polycation has no multivalent group name")
    group = name[: -len(tail)]
    if any(ch.isdigit() for ch in group):
        group = f"({group})"
    return f"{group}{multiplying_prefix(count, compound=True)}{word}"
