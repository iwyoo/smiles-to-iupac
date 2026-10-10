"""Cations of mononuclear hydrides formed by loss of a hydride ion (P-73.2.2.1): the parent hydride name with the ending
'ylium', its substituents cited in front: azanylium, phosphanylium, phenylsulfanylium, triphenylsilylium,
chloranylium. The cationic atom has one bond fewer than its standard bonding number."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._hetero_prefixes import EXTENDED_PREFIXES
from ._substituents import format_mononuclear_prefixes, name_branch

_STEM = {
    5: ("boran", 3), 7: ("azan", 3), 8: ("oxidan", 2), 14: ("sil", 4), 15: ("phosphan", 3), 16: ("sulfan", 2),
    17: ("chloran", 1), 32: ("germ", 4), 33: ("arsan", 3), 34: ("selan", 2), 35: ("broman", 1), 50: ("stann", 4),
    51: ("stiban", 3), 52: ("tellan", 2), 53: ("iodan", 1), 82: ("plumb", 4), 83: ("bismuthan", 3),
}


def _center(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or charged[0].GetFormalCharge() not in (1, 2) or len(Chem.GetMolFrags(mol)) != 1:
        return None
    atom = charged[0]
    if atom.GetAtomicNum() not in _STEM or atom.IsInRing() or atom.GetIsotope():
        return None
    standard = _STEM[atom.GetAtomicNum()][1]
    bonds = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()
    if bonds != standard - atom.GetFormalCharge():
        return None
    if atom.GetFormalCharge() == 2 and (atom.GetAtomicNum() != 7 or atom.GetDegree() != 1 or atom.GetTotalNumHs()):
        return None
    if atom.GetAtomicNum() == 7 and atom.GetDegree() and atom.GetFormalCharge() == 1:
        return None
    if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds() if atom.GetFormalCharge() == 1) or any(
        b.GetBondTypeAsDouble() > 2.0 for b in atom.GetBonds()
    ):
        return None
    if atom.GetAtomicNum() != 8 and any(n.GetAtomicNum() == 7 for n in atom.GetNeighbors()):
        return None
    if any(n.GetAtomicNum() not in (6, 7, 8, 9, 16, 17, 34, 35, 52, 53) or n.GetAtomicNum() == atom.GetAtomicNum() for n in atom.GetNeighbors()):
        return None
    if any(a.GetIsotope() for a in mol.GetAtoms()) or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms() if a.GetIdx() != atom.GetIdx()):
        return None
    return atom


def has_hydride_ylium_shape(mol) -> bool:
    return _center(mol) is not None


def name_hydride_ylium(mol) -> str:
    center = _center(mol)
    if center is None:
        raise UnsupportedStructure("not the cation of a mononuclear hydride")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [
        name_branch(graph, n, center.GetIdx(), halogens, aromatic, mol=mol) for n in graph[center.GetIdx()]
    ]
    z = center.GetAtomicNum()
    if z == 8 and len(entries) == 1:
        return _oxylium(mol, graph, center, halogens, aromatic)
    stem = _STEM[z][0]
    parent = stem + ("ylium" if center.GetFormalCharge() == 1 else "ebis(ylium)")
    return format_mononuclear_prefixes(entries) + parent if entries else parent


def _oxylium(mol, graph, center, halogens, aromatic):
    """P-73.2.3.3: the cation of a hydroxy group that has lost its hydride is the oxy group name with 'ylium': methoxylium,
    phenoxylium, (chloroacetyl)oxylium."""
    from ._hetero_prefixes import _alkoxy, _group_names, _enclose, is_functional_carbon

    (carbon,) = graph[center.GetIdx()]
    if mol.GetAtomWithIdx(carbon).GetAtomicNum() == 7:
        from .core import smiles_to_iupac

        editable = Chem.RWMol(mol)
        editable.RemoveAtom(center.GetIdx())
        amine = editable.GetMol()
        Chem.SanitizeMol(amine)
        name = smiles_to_iupac(Chem.MolToSmiles(amine))
        if not name.endswith("amine"):
            raise UnsupportedStructure("this aminoxylium has no amine name")
        return name[:-1] + "oxylium"
    if mol.GetAtomWithIdx(carbon).GetAtomicNum() == 16:
        name, compound = name_branch(graph, carbon, center.GetIdx(), halogens, aromatic, mol=mol)
        return _enclose(name, compound) + "oxylium"
    if is_functional_carbon(mol, carbon):
        ((acyl, compound),) = _group_names(graph, mol, [carbon], center.GetIdx(), halogens, aromatic)
        return _enclose(acyl, compound) + "oxylium"
    rname, compound = name_branch(graph, carbon, center.GetIdx(), halogens, aromatic, mol=mol)
    alkoxy, _ = _alkoxy(rname, compound)
    return alkoxy[: -len("oxy")] + "oxylium"


_ONIUM_STEM = {
    7: "azanium", 8: "oxidanium", 15: "phosphanium", 16: "sulfanium", 33: "arsanium", 34: "selanium", 51: "stibanium",
    52: "telluranium", 17: "chloranium", 35: "bromanium", 53: "iodanium", 83: "bismuthanium",
}
_ONIUM_NEIGHBOURS = (6, 7, 8, 9, 16, 17, 34, 35, 52, 53)


def _onium_center(mol):
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or charged[0].GetFormalCharge() != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    atom = charged[0]
    if atom.GetAtomicNum() not in _ONIUM_STEM or atom.IsInRing() or atom.GetIsotope():
        return None
    bonds = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()) + atom.GetTotalNumHs()
    if bonds != _STEM[atom.GetAtomicNum()][1] + 1 or not atom.GetDegree():
        return None
    triple_allowed = atom.GetAtomicNum() in (7, 8)
    if any(b.GetBondTypeAsDouble() > (3.0 if triple_allowed else 2.0) for b in atom.GetBonds()):
        return None
    if atom.GetAtomicNum() in (7, 8, 16):
        if any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors()):
            return None
    elif any(n.GetAtomicNum() not in _ONIUM_NEIGHBOURS or n.GetAtomicNum() == atom.GetAtomicNum() for n in atom.GetNeighbors()):
        return None
    if atom.GetAtomicNum() == 7 and all(b.GetBondTypeAsDouble() == 1.0 for b in atom.GetBonds()):
        return None
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms() if a.GetIdx() != atom.GetIdx()):
        return None
    return atom


def has_hydride_onium_shape(mol) -> bool:
    return _onium_center(mol) is not None


def name_hydride_onium(mol) -> str:
    """P-73.1.2.1: an oxonium, sulfonium or halonium centre with ylidene and alkyl groups: ethyl(propan-2-ylidene)oxidanium,
    acetyl(methyl)chloranium."""
    from ._dipolar import _group

    center = _onium_center(mol)
    if center is None:
        raise UnsupportedStructure("not an onium cation of a mononuclear hydride")
    graph = adjacency(mol)
    if center.GetAtomicNum() == 7 and any(b.GetBondTypeAsDouble() == 3.0 for b in center.GetBonds()):
        nitrilium = _substituted_nitrilium(mol, graph, center)
        if nitrilium is not None:
            return nitrilium
    entries = []
    token = EXTENDED_PREFIXES.set(True)
    try:
        for n in graph[center.GetIdx()]:
            name, compound = _group(mol, graph, n, center.GetIdx())
            entries.append((name, compound))
    finally:
        EXTENDED_PREFIXES.reset(token)
    return format_mononuclear_prefixes(entries) + _ONIUM_STEM[center.GetAtomicNum()]


def _substituted_nitrilium(mol, graph, center):
    """P-73.1.2.1: a nitrile whose nitrogen carries a substituent takes the cationic suffix 'nitrilium' with the N prefix
    ('N-methylacetonitrilium'); formonitrile admits no substitution (P-66.5.1.2.1), so it stays an azanium."""
    from ._cited_group import subtree
    from ._dipolar import _group
    from ._hetero_prefixes import _enclose
    from .core import smiles_to_iupac

    triple = next(b for b in center.GetBonds() if b.GetBondTypeAsDouble() == 3.0)
    carbon = triple.GetOtherAtom(center).GetIdx()
    others = [n for n in graph[center.GetIdx()] if n != carbon]
    if len(others) != 1:
        return None
    (substituent,) = others
    if graph[carbon] == [center.GetIdx()]:
        return None
    editable = Chem.RWMol(mol)
    editable.GetAtomWithIdx(center.GetIdx()).SetFormalCharge(0)
    for index in sorted(subtree(graph, substituent, center.GetIdx()), reverse=True):
        editable.RemoveAtom(index)
    nitrile = editable.GetMol()
    try:
        Chem.SanitizeMol(nitrile)
        name = smiles_to_iupac(Chem.MolToSmiles(nitrile))
    except (UnsupportedStructure, Chem.rdchem.MolSanitizeException):
        return None
    if not name.endswith("nitrile") or len(nitrile.GetSubstructMatches(Chem.MolFromSmarts("C#N"))) > 1:
        return None
    group, compound = _group(mol, graph, substituent, center.GetIdx())
    return f"N-{_enclose(group, compound)}{name[:-1]}ium"


_POLY_STEM = {8: "oxylium", 16: "sulfanylium", 34: "selanylium", 52: "tellanylium"}


def _poly_centers(mol):
    centers = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(centers) < 2 or len(Chem.GetMolFrags(mol)) != 1 or len({a.GetAtomicNum() for a in centers}) != 1:
        return None
    if any(
        a.GetFormalCharge() != 1
        or a.GetAtomicNum() not in _POLY_STEM
        or a.GetDegree() != 1
        or a.GetTotalNumHs()
        or a.GetIsotope()
        or a.GetNeighbors()[0].GetAtomicNum() != 6
        or mol.GetBondBetweenAtoms(a.GetIdx(), a.GetNeighbors()[0].GetIdx()).GetBondTypeAsDouble() != 1.0
        for a in centers
    ):
        return None
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        return None
    centre_indices = {c.GetIdx() for c in centers}
    if any(a.GetAtomicNum() not in (1, 6, 7, 8, 9, 17, 35, 53) and a.GetIdx() not in centre_indices for a in mol.GetAtoms()):
        return None
    return centers


def has_poly_ylium_shape(mol) -> bool:
    return _poly_centers(mol) is not None


def name_poly_ylium(mol) -> str:
    """P-73.5.1.2: ylium centres on identical hydroxy-type groups of one skeleton: (ethane-1,2-diyl)bis(oxylium),
    (pyridine-2,6-diyl)bis(sulfanylium). The skeleton's multivalent group name is read from the diacetate of its diol."""
    from ._numerals import multiplying_prefix
    from .core import smiles_to_iupac

    centers = _poly_centers(mol)
    if centers is None:
        raise UnsupportedStructure("not a polycation of identical ylium centres")
    editable = Chem.RWMol(mol)
    for center in centers:
        target = editable.GetAtomWithIdx(center.GetIdx())
        target.SetFormalCharge(0)
        target.SetNoImplicit(True)
        target.SetNumExplicitHs(0)
        target.SetNumRadicalElectrons(0)
        target.SetAtomicNum(8)
        acyl = editable.AddAtom(Chem.Atom(6))
        oxo = editable.AddAtom(Chem.Atom(8))
        methyl = editable.AddAtom(Chem.Atom(6))
        editable.AddBond(center.GetIdx(), acyl, Chem.BondType.SINGLE)
        editable.AddBond(acyl, oxo, Chem.BondType.DOUBLE)
        editable.AddBond(acyl, methyl, Chem.BondType.SINGLE)
    surrogate = editable.GetMol()
    Chem.SanitizeMol(surrogate)
    name = smiles_to_iupac(Chem.MolToSmiles(surrogate))
    count = len(centers)
    tail = " " + multiplying_prefix(count) + "acetate"
    if not name.endswith(tail):
        raise UnsupportedStructure("the skeleton of this polycation has no multivalent group name")
    group = name[: -len(tail)]
    if any(ch.isdigit() for ch in group):
        group = f"({group})"
    return f"{group}{multiplying_prefix(count, compound=True)}({_POLY_STEM[centers[0].GetAtomicNum()]})"
