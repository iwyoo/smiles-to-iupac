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


_DERIVATIVES = {
    "amide": Chem.MolFromSmarts("[CX3;H0](=O)[NX3;H2]"),
    "hydrazide": Chem.MolFromSmarts("[CX3;H0](=O)[NX3;H1][NX3;H2]"),
    "nitrile": Chem.MolFromSmarts("[CX2]#[NX1]"),
    "ester": Chem.MolFromSmarts("[CX3;H0](=O)[OX2;H0][CX4]"),
}
_ENDINGS = {"amide": "amide", "hydrazide": "ohydrazide", "nitrile": "onitrile"}


def _as_acid(mol):
    """(acid mol, derivative kind, alcohol group of an ester) for a polyhydroxy chain ending in an acid derivative."""
    found = [(kind, match) for kind, pattern in _DERIVATIVES.items() for match in mol.GetSubstructMatches(pattern)]
    found = [item for item in found if not (item[0] == "amide" and any(m[0] == item[1][0] for k, m in found if k == "hydrazide"))]
    if len(found) != 1:
        return None
    kind, match = found[0]
    editable = Chem.RWMol(mol)
    carbon = match[0]
    removal, ester_group = [], None
    if kind == "amide":
        editable.GetAtomWithIdx(match[2]).SetAtomicNum(8)
    elif kind == "hydrazide":
        removal.append(match[3])
        editable.GetAtomWithIdx(match[2]).SetAtomicNum(8)
    elif kind == "nitrile":
        nitrogen = match[1]
        editable.GetAtomWithIdx(nitrogen).SetAtomicNum(8)
        editable.GetBondBetweenAtoms(carbon, nitrogen).SetBondType(Chem.BondType.DOUBLE)
        hydroxy = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(carbon, hydroxy, Chem.BondType.SINGLE)
    else:
        from ._cited_group import cited_group, subtree
        from ._common import adjacency

        graph = adjacency(mol)
        try:
            ester_group = cited_group(mol, graph, match[3], match[2])[0]
        except Exception:
            return None
        removal.extend(subtree(graph, match[3], match[2]))
    for index in sorted(removal, reverse=True):
        editable.RemoveAtom(index)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    acid = editable.GetMol()
    try:
        Chem.SanitizeMol(acid)
    except Exception:
        return None
    return acid, kind, ester_group


def sugar_acid_derivative_name(mol):
    prepared = _as_acid(mol)
    if prepared is None:
        return None
    acid, kind, alcohol = prepared
    name = sugar_alcohol_acid_name(acid)
    if name is None:
        return None
    for ending in ("onic acid", "uronic acid"):
        if name.endswith(ending):
            stem = name[: -len(" acid")]
            break
    else:
        return None
    if kind == "ester":
        return f"{alcohol} {stem[:-2]}ate"
    return stem[:-2] + _ENDINGS[kind]


_RING_ESTER = Chem.MolFromSmarts("[CX3;R;H0](=O)[OX2;R;H0][CX4;R]")


def sugar_lactone_name(mol):
    """'D-glucono-1,5-lactone' for the internal ester of an aldonic acid (P-102.5.6.6.2.2): the ring is opened at the
    ester oxygen, the acid named, and the locants of the carboxy carbon and of the carbon that closes the ring cited."""
    from ._sugar_substituted import _opened_ether

    matches = mol.GetSubstructMatches(_RING_ESTER)
    if len(matches) != 1 or mol.GetRingInfo().NumRings() != 1:
        return None
    carbonyl, _, oxygen, closing = matches[0]
    try:
        opened = _opened_ether(mol, oxygen, closing)
    except (ValueError, RuntimeError):
        return None
    acid = sugar_alcohol_acid_name(opened)
    if acid is None or not acid.endswith("onic acid"):
        return None
    position, seen, frontier = {carbonyl: 1}, {carbonyl}, [carbonyl]
    while frontier:
        atom = frontier.pop(0)
        for n in opened.GetAtomWithIdx(atom).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIdx() not in seen:
                seen.add(n.GetIdx())
                position[n.GetIdx()] = position[atom] + 1
                frontier.append(n.GetIdx())
    return f"{acid[: -len('ic acid')]}o-1,{position[closing]}-lactone"


def has_sugar_lactone_shape(mol) -> bool:
    return sugar_lactone_name(mol) is not None


def name_sugar_lactone(mol) -> str:
    return sugar_lactone_name(mol)
