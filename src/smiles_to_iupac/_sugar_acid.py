"""Alditols, aldonic, uronic and aldaric acids (P-102.5.6.5, P-102.5.6.6).

The unbranched polyhydroxy chain is brought back to the aldose it derives from (a terminal CH2OH or COOH becomes CHO,
the other end CH2OH), named as that aldose, and the 'ose' ending becomes 'itol', 'onic acid', 'uronic acid' or
'aric acid'. A chain derivable from two aldoses takes the parent chosen by P-102.4 (c): the first stem alphabetically,
then D before L; a chain that is its own mirror image is a meso form and carries no D or L.
"""

from rdkit import Chem

from ._carbohydrate import has_open_chain_aldose_shape, name_open_chain_aldose

_MIN_CARBONS, _MAX_CARBONS = 4, 7


def _end_kind(mol, carbon):
    atom = mol.GetAtomWithIdx(carbon)
    oxygens = [(n, mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble()) for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
    singles = [n for n, order in oxygens if order == 1.0]
    doubles = [n for n, order in oxygens if order == 2.0]
    if len(singles) == 1 and not doubles and atom.GetTotalNumHs() == 2:
        return "CH2OH", singles[0].GetIdx(), None
    if len(doubles) == 1 and not singles and atom.GetTotalNumHs() == 1:
        return "CHO", None, doubles[0].GetIdx()
    if len(doubles) == 1 and len(singles) == 1 and atom.GetTotalNumHs() == 0:
        return "COOH", singles[0].GetIdx(), doubles[0].GetIdx()
    return None


def _chain(mol):
    """([C-1 .. C-n] atom indices from one end, [end kinds]) of an unbranched CH(OH) chain, else None."""
    if mol.GetRingInfo().NumRings() or any(a.GetAtomicNum() not in (6, 8) or a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    if not _MIN_CARBONS <= len(carbons) <= _MAX_CARBONS:
        return None
    neighbors = {c: [n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 6] for c in carbons}
    ends = [c for c in carbons if len(neighbors[c]) == 1]
    if len(ends) != 2 or any(len(v) > 2 for v in neighbors.values()):
        return None
    chain = [ends[0]]
    while len(chain) < len(carbons):
        following = [n for n in neighbors[chain[-1]] if n not in chain]
        if len(following) != 1:
            return None
        chain.append(following[0])
    for c in chain[1:-1]:
        atom = mol.GetAtomWithIdx(c)
        hydroxyls = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(hydroxyls) != 1 or atom.GetTotalNumHs() != 1 or mol.GetBondBetweenAtoms(c, hydroxyls[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
            return None
    kinds = [_end_kind(mol, chain[0]), _end_kind(mol, chain[-1])]
    if None in kinds or mol.GetNumAtoms() != len(chain) + sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8):
        return None
    return chain, kinds


def _as_aldose(mol, chain, kinds, reverse):
    """Name of the aldose reached by making the first (or, reversed, last) chain end the aldehyde and the other CH2OH."""
    order = list(reversed(chain)) if reverse else list(chain)
    kind_by_end = dict(zip((chain[0], chain[-1]), kinds))
    first, last = order[0], order[-1]
    editable = Chem.RWMol(mol)
    first_kind, last_kind = kind_by_end[first], kind_by_end[last]
    drop = []
    if first_kind[0] == "COOH":
        drop.append(first_kind[1])
    elif first_kind[0] == "CH2OH":
        editable.GetBondBetweenAtoms(first, first_kind[1]).SetBondType(Chem.BondType.DOUBLE)
    if last_kind[0] == "COOH":
        drop.append(last_kind[2])
    elif last_kind[0] == "CHO":
        editable.GetBondBetweenAtoms(last, last_kind[2]).SetBondType(Chem.BondType.SINGLE)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    for index in sorted(drop, reverse=True):
        editable.RemoveAtom(index)
    model = editable.GetMol()
    try:
        Chem.SanitizeMol(model)
    except Exception:
        return None
    return name_open_chain_aldose(model) if has_open_chain_aldose_shape(model) else None


def _class(kinds):
    ends = tuple(sorted(k[0] for k in kinds))
    return {
        ("CH2OH", "CH2OH"): "alditol",
        ("CH2OH", "COOH"): "aldonic",
        ("CHO", "COOH"): "uronic",
        ("COOH", "COOH"): "aldaric",
    }.get(ends)


def _suffixes(kind):
    return {"alditol": "itol", "aldonic": "onic acid", "uronic": "uronic acid", "aldaric": "aric acid"}[kind]


def sugar_alcohol_acid_name(mol):
    found = _chain(mol)
    if found is None:
        return None
    chain, kinds = found
    kind = _class(kinds)
    if kind is None:
        return None
    if kind == "aldonic":
        reverses = [kinds[0][0] != "COOH"]
    elif kind == "uronic":
        reverses = [kinds[0][0] != "CHO"]
    else:
        reverses = [False, True]
    names = [_as_aldose(mol, chain, kinds, reverse) for reverse in reverses]
    if None in names:
        return None
    if len(names) == 2 and names[0][2:] == names[1][2:] and names[0][0] != names[1][0]:
        name = names[0][2:]
    else:
        name = min(names, key=lambda n: (n[2:], n[0] != "D"))
    if not name.endswith("ose"):
        return None
    return name[: -len("ose")] + _suffixes(kind)


def has_sugar_alcohol_acid_shape(mol) -> bool:
    return sugar_alcohol_acid_name(mol) is not None


def name_sugar_alcohol_acid(mol) -> str:
    return sugar_alcohol_acid_name(mol)
