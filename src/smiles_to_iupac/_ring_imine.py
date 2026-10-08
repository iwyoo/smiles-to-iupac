"""Imines whose C=N carbon is a ring atom of a heterocycle, a partly hydrogenated ring or a fused system (P-62.3.1.1):
thiolan-2-imine, naphthalen-2(1H)-imine, N,N'-dimethylnaphthalene-1,4-diimine. The suffix 'imine' replaces 'one' on the
same skeleton with the same locants and hydro prefixes, so the =N-R groups are named as =O on a copy of the structure
and the N-substituents are cited as N-prefixes."""

import re

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import name_branch

_ONE_LOCANTS = re.compile(r"(\d+(?:,\d+)*)(?:\([0-9a-zH,]+\))?-(?:di|tri|tetra|penta|hexa)one$")
_MULTIPLIER = {1: "", 2: "di", 3: "tri", 4: "tetra"}


def _imines(mol):
    """[(ring carbon, nitrogen)] of every ring C=N; None when any other atom or bond kind is present."""
    found = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondTypeAsDouble() != 2.0 or {a.GetAtomicNum(), b.GetAtomicNum()} != {6, 7}:
            continue
        carbon, nitrogen = (a, b) if a.GetAtomicNum() == 6 else (b, a)
        if not carbon.IsInRing() or nitrogen.IsInRing():
            return None
        found.append((carbon, nitrogen))
    return found


def _ketone_of(mol, imines):
    """Copy of `mol` with each =N-R turned into =O, and the substituent roots removed from the copy."""
    editable = Chem.RWMol(mol)
    drop = set()
    for _, nitrogen in imines:
        editable.GetAtomWithIdx(nitrogen.GetIdx()).SetAtomicNum(8)
        editable.GetAtomWithIdx(nitrogen.GetIdx()).SetNumExplicitHs(0)
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetIdx() in {c.GetIdx() for c, _ in imines}:
                continue
            editable.RemoveBond(nitrogen.GetIdx(), neighbor.GetIdx())
            drop |= {i for i in _side(mol, neighbor.GetIdx(), nitrogen.GetIdx())}
    for idx in sorted(drop, reverse=True):
        editable.RemoveAtom(idx)
    ketone = editable.GetMol()
    Chem.SanitizeMol(ketone)
    return ketone


def _side(mol, start, blocked):
    seen, stack = set(), [start]
    while stack:
        idx = stack.pop()
        if idx in seen:
            continue
        seen.add(idx)
        stack.extend(n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetIdx() != blocked)
    return seen


def _name(mol):
    from .core import smiles_to_iupac

    imines = _imines(mol)
    if not imines or len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("no ring imine")
    nitrogens = {n.GetIdx() for _, n in imines}
    if any(
        a.GetIdx() not in nitrogens
        and (a.GetFormalCharge() or a.GetIsotope() or (not a.IsInRing() and a.GetAtomicNum() not in (6, *HALOGEN_PREFIXES)))
        for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("charged or isotopic atoms and acyclic atoms other than carbon, halogens and the imino nitrogens are not supported yet")
    if any(n.GetFormalCharge() or n.GetDegree() > 2 for _, n in imines):
        raise UnsupportedStructure("an imino nitrogen with more than one substituent is not supported")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    cited = []
    for carbon, nitrogen in imines:
        roots = [r for r in graph[nitrogen.GetIdx()] if r != carbon.GetIdx()]
        if any(mol.GetBondBetweenAtoms(nitrogen.GetIdx(), r).GetBondTypeAsDouble() != 1.0 for r in roots):
            raise UnsupportedStructure("a multiple bond on an imino nitrogen is not supported yet")
        if roots:
            substituent, compound = name_branch(graph, roots[0], nitrogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
            cited.append(f"({substituent})" if compound else substituent)
        else:
            cited.append(None)
    if len(set(cited)) != 1:
        raise UnsupportedStructure("imino groups with different N-substituents are not supported yet")
    ketone_name = smiles_to_iupac(Chem.MolToSmiles(_ketone_of(mol, imines)))
    if not ketone_name.endswith("one") or ketone_name.endswith(("thione", "selenone", "tellurone")):
        raise UnsupportedStructure("the ketone of this ring imine is not named with the suffix 'one'")
    name = ketone_name[:-3] + "imine"
    if cited[0] is None:
        return name
    if not ketone_name[0].isalpha():
        raise UnsupportedStructure("N-substituents beside other prefixes on a ring imine are not supported yet")
    locants = _ONE_LOCANTS.search(ketone_name)
    count = len(imines)
    prefix_locants = ",".join(f"N{loc}" for loc in locants.group(1).split(",")) if locants else "N"
    word = cited[0] if count == 1 else _MULTIPLIER[count] + cited[0]
    return f"{prefix_locants}-{word}{name}"


def has_ring_imine_shape(mol) -> bool:
    try:
        _name(mol)
    except Exception:
        return False
    return True


def name_ring_imine(mol) -> str:
    return _name(mol)
