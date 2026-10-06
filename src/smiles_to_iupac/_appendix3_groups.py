"""Characteristic groups on an Appendix 3 parent (P-41, P-44.1.1, P-65, P-66, P-101.7.1): which class is principal,
where its members sit on the skeleton and which neighbours are cited as prefixes through `name_branch`."""

from collections import Counter
from dataclasses import dataclass, field

from ._common import UnsupportedStructure

SENIORITY = (
    "acid",
    "sulfonic",
    "anhydride",
    "ester",
    "halide",
    "amide",
    "sulfonamide",
    "amidine",
    "hydrazide",
    "nitrile",
    "aldehyde",
    "ketone",
    "alcohol",
    "thiol",
    "amine",
    "imine",
)
_HALOGENS = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_OTHER_ELEMENTS = frozenset({1, 6, 7, 8, 9, 16, 17, 35, 53})


@dataclass
class Group:
    cls: str
    label: str
    kind: str = ""
    atoms: tuple = ()
    root: int = -1
    extra: dict = field(default_factory=dict)


def _double_oxygens(mol, carbon, excluded=()):
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(carbon).GetNeighbors()
        if n.GetIdx() not in excluded
        and n.GetAtomicNum() == 8
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]


def _single_neighbors(mol, carbon, excluded):
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(carbon).GetNeighbors()
        if n.GetIdx() not in excluded and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0
    ]


def acyl_class(mol, carbon, anchor, mapped):
    """(class, detail) for the carbon of a C(=O)X / C#N / CH=O group whose other side is `anchor`, else None."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() and carbon in mapped:
        return None
    triple = [
        n.GetIdx()
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 7 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 3.0
    ]
    if triple and atom.GetDegree() == 2:
        return "nitrile", {"atoms": (carbon, triple[0])}
    oxo = _double_oxygens(mol, carbon)
    if not oxo:
        return _amidine(mol, carbon, anchor, mapped)
    if len(oxo) != 1:
        return None
    rest = _single_neighbors(mol, carbon, {anchor, oxo[0]})
    if not rest:
        if atom.GetDegree() == 2 and atom.GetTotalNumHs() == 1:
            return "aldehyde", {"atoms": (carbon, oxo[0])}
        return None
    if len(rest) != 1 or atom.GetDegree() != 3:
        return None
    x = mol.GetAtomWithIdx(rest[0])
    if x.GetAtomicNum() == 8 and x.GetDegree() == 1 and x.GetTotalNumHs() == 1:
        return "acid", {"atoms": (carbon, oxo[0], rest[0])}
    if x.GetAtomicNum() == 8 and x.GetDegree() == 2:
        alkyl = next(n for n in x.GetNeighbors() if n.GetIdx() != carbon)
        if alkyl.GetIdx() in mapped or alkyl.GetAtomicNum() not in (6,):
            raise UnsupportedStructure("a lactone or non-carbon ester group on a natural-product skeleton is not supported")
        if _double_oxygens(mol, alkyl.GetIdx()):
            return "anhydride", {"atoms": (carbon, oxo[0], rest[0]), "oxygen": rest[0], "other": alkyl.GetIdx()}
        return "ester", {"atoms": (carbon, oxo[0], rest[0]), "oxygen": rest[0], "alkyl": alkyl.GetIdx()}
    if x.GetAtomicNum() in _HALOGENS and x.GetDegree() == 1:
        return "halide", {"atoms": (carbon, oxo[0], rest[0]), "halide": _HALOGENS[x.GetAtomicNum()]}
    if x.GetAtomicNum() == 7 and not x.IsInRing() and x.GetFormalCharge() == 0:
        hydrazine = [n for n in x.GetNeighbors() if n.GetAtomicNum() == 7]
        if hydrazine:
            terminal = hydrazine[0]
            if len(hydrazine) > 1 or terminal.IsInRing() or terminal.GetFormalCharge() or terminal.GetIsAromatic():
                return None
            return "hydrazide", {"atoms": (carbon, oxo[0], rest[0], terminal.GetIdx()), "nitrogen": rest[0], "terminal": terminal.GetIdx()}
        return "amide", {"atoms": (carbon, oxo[0], rest[0]), "nitrogen": rest[0]}
    return None


def _amidine(mol, carbon, anchor, mapped):
    """R-C(=NR')-NR2 on a carbon that carries no oxygen (P-66.4.1.1)."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetDegree() != 3:
        return None
    imino = amino = None
    for n in atom.GetNeighbors():
        if n.GetIdx() == anchor or n.GetIdx() in mapped or n.GetAtomicNum() != 7 or n.IsInRing() or n.GetIsAromatic() or n.GetFormalCharge():
            continue
        order = mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble()
        if order == 2.0 and imino is None:
            imino = n.GetIdx()
        elif order == 1.0 and amino is None:
            amino = n.GetIdx()
    if imino is None or amino is None:
        return None
    return "amidine", {"atoms": (carbon, imino, amino), "nitrogen": amino, "imino": imino}


def _o_acyl(mol, oxygen, parent, mapped):
    """The acyl carbon of an R-C(=O)-O- group on `parent`, else None."""
    atom = mol.GetAtomWithIdx(oxygen)
    if atom.GetAtomicNum() != 8 or atom.GetDegree() != 2:
        return None
    acyl = next(n for n in atom.GetNeighbors() if n.GetIdx() != parent)
    if acyl.GetAtomicNum() != 6 or acyl.GetIdx() in mapped:
        return None
    if len(_double_oxygens(mol, acyl.GetIdx())) != 1:
        return None
    return acyl.GetIdx()


def _amine_nitrogen(mol, nitrogen, anchor, mapped):
    atom = mol.GetAtomWithIdx(nitrogen)
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() or atom.GetIsAromatic() or atom.IsInRing():
        return None
    if any(mol.GetBondBetweenAtoms(nitrogen, n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in atom.GetNeighbors()):
        return None
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != anchor]
    for n in others:
        if n.GetIdx() in mapped:
            return None
        if n.GetAtomicNum() != 6:
            return None
        if n.GetAtomicNum() == 6 and (_double_oxygens(mol, n.GetIdx()) or any(
            m.GetAtomicNum() == 16 and m.GetDegree() > 2 for m in n.GetNeighbors()
        )):
            return None
    return [n.GetIdx() for n in others]


def classify(mol, mapping, terminals):
    """(groups, branches, attach) for the atoms outside the skeleton `mapping` ({label: atom}); `attach` lists the
    (label, dummy atom) of each free valence, or is None.

    `terminals`: labels of skeleton chain atoms that end a chain (a functional carbon there is part of the parent)."""
    mapped = set(mapping.values())
    groups, branches, attach = [], [], []
    for label, atom_idx in mapping.items():
        atom = mol.GetAtomWithIdx(atom_idx)
        outside = [n for n in atom.GetNeighbors() if n.GetIdx() not in mapped]
        dummy = [n.GetIdx() for n in outside if n.GetAtomicNum() == 0]
        if dummy:
            if any(mol.GetBondBetweenAtoms(atom_idx, d).GetBondTypeAsDouble() != 1.0 for d in dummy):
                raise UnsupportedStructure("a multiple-bond free valence on an Appendix 3 parent is not supported")
            attach.extend((label, d) for d in dummy)
            outside = [n for n in outside if n.GetAtomicNum() != 0]
        if label.startswith("_"):
            if outside:
                raise UnsupportedStructure("a substituent on an unnumbered skeleton atom is not supported")
            continue
        if label in terminals and not atom.IsInRing():
            inside = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in mapped]
            found = acyl_class(mol, atom_idx, inside[0], mapped) if outside and len(inside) == 1 else None
            if found is not None:
                cls, detail = found
                groups.append(Group(cls, label, "o", detail["atoms"], -1, detail))
                continue
        for neighbor in outside:
            idx = neighbor.GetIdx()
            order = mol.GetBondBetweenAtoms(atom_idx, idx).GetBondTypeAsDouble()
            number = neighbor.GetAtomicNum()
            if number == 8 and neighbor.GetDegree() == 1 and order == 2.0:
                groups.append(Group("ketone", label, "", (idx,), idx))
            elif number == 8 and neighbor.GetDegree() == 1 and order == 1.0 and neighbor.GetTotalNumHs() == 1:
                groups.append(Group("alcohol", label, "", (idx,), idx))
            elif number == 16 and neighbor.GetDegree() == 1 and order == 1.0 and neighbor.GetTotalNumHs() == 1:
                groups.append(Group("thiol", label, "", (idx,), idx))
            elif number == 7 and order == 1.0 and (substituents := _amine_nitrogen(mol, idx, atom_idx, mapped)) is not None:
                groups.append(Group("amine", label, "", (idx,), idx, {"substituents": substituents}))
            elif number == 8 and order == 1.0 and _o_acyl(mol, idx, atom_idx, mapped) is not None:
                acyl = _o_acyl(mol, idx, atom_idx, mapped)
                groups.append(Group("ester_o", label, "", (idx,), idx, {"acyl": acyl}))
            elif number == 6 and order == 1.0 and (found := acyl_class(mol, idx, atom_idx, mapped)) is not None:
                cls, detail = found
                groups.append(Group(cls, label, "c", detail["atoms"], idx, detail))
            elif number == 16 and order == 1.0 and _sulfonic(mol, idx, atom_idx):
                groups.append(Group("sulfonic", label, "c", (idx,), idx))
            elif number == 16 and order == 1.0 and (nitrogen := _sulfonamide_nitrogen(mol, idx, mapped)) is not None:
                groups.append(Group("sulfonamide", label, "c", (idx,), idx, {"nitrogen": nitrogen}))
            elif number == 7 and order == 2.0 and (substituents := _imine_substituents(mol, idx, atom_idx, mapped)) is not None:
                groups.append(Group("imine", label, "", (idx,), idx, {"substituents": substituents}))
            else:
                branches.append((label, idx))
    return groups, branches, tuple(attach) or None


def _sulfonic(mol, sulfur, parent):
    atom = mol.GetAtomWithIdx(sulfur)
    if atom.GetDegree() != 4:
        return False
    oxo = _double_oxygens(mol, sulfur)
    hydroxy = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1]
    return len(oxo) == 2 and len(hydroxy) == 1


def _sulfonamide_nitrogen(mol, sulfur, mapped):
    atom = mol.GetAtomWithIdx(sulfur)
    if atom.GetDegree() != 4 or len(_double_oxygens(mol, sulfur)) != 2:
        return None
    nitrogens = [
        n for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 7 and not n.IsInRing() and not n.GetFormalCharge() and n.GetIdx() not in mapped
    ]
    return nitrogens[0].GetIdx() if len(nitrogens) == 1 else None


def _imine_substituents(mol, nitrogen, anchor, mapped):
    """Atoms on the nitrogen of a C=N-R group (R may be H); None when R is not a carbon or oxygen group."""
    atom = mol.GetAtomWithIdx(nitrogen)
    if atom.IsInRing() or atom.GetIsAromatic() or atom.GetFormalCharge():
        return None
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != anchor]
    for n in others:
        if n.GetIdx() in mapped or n.GetAtomicNum() not in (6, 8):
            return None
        if mol.GetBondBetweenAtoms(nitrogen, n.GetIdx()).GetBondTypeAsDouble() != 1.0:
            return None
    return [n.GetIdx() for n in others]


def branch_counts(mol, root, anchor, mapped):
    """Counter of principal-capable group classes inside the side branch rooted at `root` (reached from `anchor`)."""
    seen, stack, counts = {anchor}, [root], Counter()
    while stack:
        idx = stack.pop()
        if idx in seen or idx in mapped:
            continue
        seen.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        number = atom.GetAtomicNum()
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        if number == 6:
            found = acyl_class(mol, idx, next((n for n in neighbors if n in seen), -1), mapped)
            if found is not None:
                counts[found[0]] += 1
            elif _double_oxygens(mol, idx) and atom.GetDegree() == 3:
                counts["ketone"] += 1
        elif number == 8 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1:
            counts["alcohol"] += 1
        elif number == 16 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1:
            counts["thiol"] += 1
        elif number == 7 and not atom.GetIsAromatic() and not atom.IsInRing() and not atom.GetFormalCharge():
            carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            single = all(mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 1.0 for n in atom.GetNeighbors())
            if len(carbons) == atom.GetDegree() and single and not any(_double_oxygens(mol, c.GetIdx()) for c in carbons):
                counts["amine"] += 1
        stack.extend(n for n in neighbors if n not in seen)
    return counts


def reject_exotic(mol, mapped):
    for atom in mol.GetAtoms():
        if atom.GetIdx() in mapped:
            continue
        number = atom.GetAtomicNum()
        if number == 0:
            continue
        if number not in _OTHER_ELEMENTS:
            raise UnsupportedStructure("an element outside C, H, N, O, S and the halogens is not supported here")
        if number == 16 and atom.GetDegree() > 2 and not _sulfonic(mol, atom.GetIdx(), -1) and _sulfonamide_nitrogen(mol, atom.GetIdx(), mapped) is None:
            raise UnsupportedStructure("a sulfur group other than sulfanyl or sulfo is not supported here")
        if atom.GetFormalCharge() and not (number == 7 and sum(n.GetAtomicNum() == 8 for n in atom.GetNeighbors()) == 2):
            raise UnsupportedStructure("a charged group is not supported here")
