"""Alditols, aldonic, uronic and aldaric acids (P-102.5.6.5, P-102.5.6.6).

The unbranched polyhydroxy chain is brought back to the aldose it derives from (a terminal CH2OH or COOH becomes CHO,
the other end CH2OH), named as that aldose, and the 'ose' ending becomes 'itol', 'onic acid', 'uronic acid' or
'aric acid'. A chain derivable from two aldoses takes the parent chosen by P-102.4 (c): the first stem alphabetically,
then D before L; a chain that is its own mirror image is a meso form and carries no D or L.
"""

import re

from rdkit import Chem

from ._carbohydrate import has_open_chain_aldose_shape, name_open_chain_aldose
from ._common import UnsupportedStructure, adjacency, alphanumerical_name_key

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


def _resolved(mol):
    """(class, name, chains) of an alditol, aldonic, uronic or aldaric acid: each chain lists the carbons from C-1 for an
    orientation that gives the chosen name (both of them when the name is the same either way)."""
    found = _chain(mol)
    if found is None:
        return None
    chain, kinds = found
    kind = _class(kinds)
    if kind is None or (kind == "aldaric" and len(chain) == 4):
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
        name, winners = names[0][2:], reverses
    else:
        name = min(names, key=lambda n: (n[2:], n[0] != "D"))
        winners = [reverse for reverse, n in zip(reverses, names) if n == name]
    if not name.endswith("ose"):
        return None
    return kind, name[: -len("ose")] + _suffixes(kind), [chain[::-1] if reverse else chain for reverse in winners]


def sugar_alcohol_acid_name(mol):
    resolved = _resolved(mol)
    return resolved[1] if resolved else None


def has_sugar_alcohol_acid_shape(mol) -> bool:
    return sugar_alcohol_acid_name(mol) is not None


def name_sugar_alcohol_acid(mol) -> str:
    return sugar_alcohol_acid_name(mol)


_ACID_CARBON = "[CX3;H0;$(C-[CX4;$(C-[#8,#7,#9,#17,#35,#53])])]"
_DERIVATIVES = {
    "amide": Chem.MolFromSmarts(f"{_ACID_CARBON}(=O)[NX3;H2]"),
    "hydrazide": Chem.MolFromSmarts(f"{_ACID_CARBON}(=O)[NX3;H1][NX3;H2]"),
    "nitrile": Chem.MolFromSmarts("[CX2;$(C-[CX4;$(C-[#8,#7,#9,#17,#35,#53])])]#[NX1]"),
    "ester": Chem.MolFromSmarts(f"{_ACID_CARBON}(=O)[OX2;H0][CX4]"),
    "halide": Chem.MolFromSmarts(f"{_ACID_CARBON}(=O)[F,Cl,Br,I;X1]"),
}
_HALIDE_WORDS = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
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
    elif kind == "halide":
        ester_group = _HALIDE_WORDS[mol.GetAtomWithIdx(match[2]).GetAtomicNum()]
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
    for index in sorted(set(removal), reverse=True):
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


_OXO_SUBSTITUENT = {
    "amide": "{n}-amino-{n}-deoxy-{n}-oxo-",
    "hydrazide": "{n}-deoxy-{n}-hydrazinyl-{n}-oxo-",
}


def _aldaric_with_modified_ends(mol):
    """(aldaric acid mol, {end carbon: (kind, ester group)}) for a chain whose carboxy ends are esters or amides, else None."""
    from ._cited_group import cited_group, subtree

    if mol.GetRingInfo().NumRings():
        return None
    ends = [(kind, match) for kind in ("ester", "hydrazide", "amide") for match in mol.GetSubstructMatches(_DERIVATIVES[kind])]
    ends = [(kind, match) for kind, match in ends if kind != "amide" or all(m[0] != match[0] for k, m in ends if k == "hydrazide")]
    if not 1 <= len(ends) <= 2:
        return None
    graph = adjacency(mol)
    editable = Chem.RWMol(mol)
    removal, modified = [], {}
    for kind, match in ends:
        if kind in ("amide", "hydrazide"):
            editable.GetAtomWithIdx(match[2]).SetAtomicNum(8)
            if kind == "hydrazide":
                removal.append(match[3])
            modified[match[0]] = (kind, None)
            continue
        try:
            group = cited_group(mol, graph, match[3], match[2])[0]
        except Exception:
            return None
        removal.extend(subtree(graph, match[3], match[2]))
        modified[match[0]] = (kind, group)
    removed = set(removal)
    if len(removed) != len(removal) or removed & {match[0] for _, match in ends}:
        return None
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    acid = editable.GetMol()
    try:
        Chem.SanitizeMol(acid)
    except Exception:
        return None
    kept = [i for i in range(mol.GetNumAtoms()) if i not in removed]
    new_index = {old: new for new, old in enumerate(kept)}
    return acid, {new_index[carbon]: entry for carbon, entry in modified.items()}


def _aldonic_of_amic_end(acid, amide_carbon):
    """The aldonic acid reached by reducing the carboxy group at `amide_carbon` to CH2OH."""
    editable = Chem.RWMol(acid)
    oxygen = next(
        n.GetIdx()
        for n in acid.GetAtomWithIdx(amide_carbon).GetNeighbors()
        if acid.GetBondBetweenAtoms(amide_carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    )
    editable.RemoveAtom(oxygen)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    aldonic = editable.GetMol()
    Chem.SanitizeMol(aldonic)
    return sugar_alcohol_acid_name(aldonic)


def aldaric_derivative_name(mol):
    """Half esters, amic acids, ester amides and diamides of an aldaric acid (P-102.5.6.6.5.3): '6-methyl hydrogen
    D-glucarate', '6-amino-6-deoxy-6-oxo-D-gluconic acid', 'methyl 6-amino-6-deoxy-6-oxo-D-gluconate', 'D-glucaramide'."""
    prepared = _aldaric_with_modified_ends(mol)
    if prepared is None:
        return None
    acid, modified = prepared
    resolved = _resolved(acid)
    if resolved is None or resolved[0] != "aldaric":
        return None
    _, name, chains = resolved
    first, last = chains[0][0], chains[0][-1]
    kinds = {end: modified.get(end, ("acid", None)) for end in (first, last)}
    stem = name[: -len("ic acid")]
    size = len(chains[0])
    by_kind = sorted(kind for kind, _ in kinds.values())
    if by_kind == ["amide", "amide"]:
        return stem + "amide"
    if by_kind == ["hydrazide", "hydrazide"]:
        return stem + "odihydrazide"
    if by_kind == ["ester", "ester"]:
        groups = [kinds[first][1], kinds[last][1]]
        if groups[0] == groups[1]:
            group = f"({groups[0]})" if re.search(r"[\d\-\s()\[\]]", groups[0]) else groups[0]
            return f"di{group} {stem}ate"
        cited = min(
            sorted(((kinds[chain[0]][1], 1), (kinds[chain[-1]][1], size)), key=lambda part: alphanumerical_name_key(part[0]))
            for chain in chains
        )
        return " ".join(f"{locant}-{group}" for group, locant in cited) + f" {stem}ate"
    if by_kind == ["acid", "ester"]:
        (ester_end,) = (end for end, (kind, _) in kinds.items() if kind == "ester")
        locant = min(1 if chain[0] == ester_end else size for chain in chains)
        return f"{locant}-{kinds[ester_end][1]} hydrogen {stem}ate"
    senior_amide = by_kind == ["amide", "hydrazide"]
    junior_kind = "hydrazide" if senior_amide else next((kind for kind in _OXO_SUBSTITUENT if kind in by_kind), None)
    if junior_kind is None or not (senior_amide or {"acid", "ester"} & set(by_kind)):
        return None
    junior = next(end for end, (kind, _) in kinds.items() if kind == junior_kind)
    aldonic = _aldonic_of_amic_end(acid, junior)
    if aldonic is None or not aldonic.endswith("onic acid"):
        return None
    prefix = _OXO_SUBSTITUENT[junior_kind].format(n=size)
    if senior_amide:
        return f"{prefix}{aldonic[: -len('ic acid')]}amide"
    if "acid" in by_kind:
        return prefix + aldonic
    (ester_end,) = (end for end, (kind, _) in kinds.items() if kind == "ester")
    return f"{kinds[ester_end][1]} {prefix}{aldonic[: -len('ic acid')]}ate"


def sugar_acid_derivative_name(mol):
    aldaric = aldaric_derivative_name(mol)
    if aldaric is not None:
        return aldaric
    prepared = _as_acid(mol)
    if prepared is None:
        return None
    acid, kind, alcohol = prepared
    name = sugar_alcohol_acid_name(acid) or substituted_chain_name(acid)
    if name is None:
        return None
    for ending in ("onic acid", "uronic acid"):
        if name.endswith(ending):
            stem = name[: -len(" acid")]
            break
    else:
        return None
    if kind == "halide":
        return f"{stem[:-2]}oyl {alcohol}"
    if kind == "ester":
        return f"{alcohol} {stem[:-2]}ate"
    return stem[:-2] + _ENDINGS[kind]


_CARBOXYLATE = Chem.MolFromSmarts("[CX3;H0](=O)[OX1;-1]")


def sugar_carboxylate_name(mol):
    """P-72.2.2.2.1.1, P-102.5.6.6.2.1, P-102.5.6.6.4: the carboxylate of an aldonic, uronic or glycosiduronic acid keeps
    the sugar-style name of the acid, whose 'ic acid' ending becomes 'ate' ('D-gluconate', 'α-D-glucopyranuronate')."""
    from ._glycoside import has_glycoside_shape, name_glycoside
    from .core import smiles_to_iupac

    anions = mol.GetSubstructMatches(_CARBOXYLATE)
    if len(anions) != 1 or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() and a.GetIdx() != anions[0][2] for a in mol.GetAtoms()):
        return None
    neutral = Chem.RWMol(mol)
    oxygen = neutral.GetAtomWithIdx(anions[0][2])
    oxygen.SetFormalCharge(0)
    oxygen.SetNoImplicit(False)
    acid = neutral.GetMol()
    Chem.SanitizeMol(acid)
    for has_shape, name_function in (
        (has_sugar_ring_acid_shape, name_sugar_ring_acid),
        (has_ketoaldonic_chain_shape, name_ketoaldonic_chain),
        (has_glycoside_shape, name_glycoside),
        (has_sugar_alcohol_acid_shape, name_sugar_alcohol_acid),
        (has_substituted_chain_shape, name_substituted_chain),
    ):
        if has_shape(acid):
            name = name_function(acid)
            return name[: -len("ic acid")] + "ate" if name.endswith(("onic acid", "uronic acid")) else None
    return None


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


_CHAIN_HALOGENS = {9, 17, 35, 53}
_RETAINED_ALDITOLS = {"L-galactose": "L-fucitol", "L-mannose": "L-rhamnitol"}


def _decorated_chain(mol, allow_oxo=False):
    """([C-1 .. C-n] from one end, per-carbon (kind, exo atom, root)) of an unbranched carbon chain that carries
    hydroxy groups, deoxy positions, ethers, amino groups or halogens, else None."""
    from ._sugar_substituted import _CHALCOGEN_WORDS, _HALOGENS  # noqa: F401

    if mol.GetRingInfo().NumRings() or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    carbon_graph = {a.GetIdx(): [n for n in graph[a.GetIdx()] if mol.GetAtomWithIdx(n).GetAtomicNum() == 6] for a in mol.GetAtoms() if a.GetAtomicNum() == 6}
    pieces = []
    seen = set()
    for start in carbon_graph:
        if start in seen:
            continue
        stack, piece = [start], set()
        while stack:
            current = stack.pop()
            if current in piece:
                continue
            piece.add(current)
            stack.extend(carbon_graph[current])
        seen |= piece
        pieces.append(piece)
    candidates = [p for p in pieces if _MIN_CARBONS <= len(p) <= _MAX_CARBONS]
    if not candidates:
        return None

    def oxygenation(piece):
        return sum(1 for c in piece for n in graph[c] if n not in piece and mol.GetAtomWithIdx(n).GetAtomicNum() in (7, 8))

    best = max(oxygenation(p) for p in candidates)
    main = [p for p in candidates if oxygenation(p) == best]
    if len(main) != 1:
        return None
    piece = main[0]
    ends = [c for c in piece if len(carbon_graph[c]) == 1]
    if len(ends) != 2 or any(len(carbon_graph[c]) > 2 for c in piece):
        return None
    chain = [ends[0]]
    while len(chain) < len(piece):
        following = [n for n in carbon_graph[chain[-1]] if n not in chain]
        if len(following) != 1:
            return None
        chain.append(following[0])
    if any(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0 for a, b in zip(chain, chain[1:])):
        return None
    decorations = []
    for carbon in chain:
        exo = [n for n in graph[carbon] if n not in piece]
        if len(exo) == 2 and carbon in (ends[0], ends[-1]) and all(
            mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1 for n in exo
        ) and sorted(mol.GetBondBetweenAtoms(carbon, n).GetBondTypeAsDouble() for n in exo) == [1.0, 2.0]:
            decorations.append(("COOH", None, None))
            continue
        if len(exo) > 1:
            return None
        if not exo:
            decorations.append(("H", None, None))
            continue
        (x,) = exo
        atom = mol.GetAtomWithIdx(x)
        z = atom.GetAtomicNum()
        if allow_oxo and z == 8 and atom.GetDegree() == 1 and mol.GetBondBetweenAtoms(carbon, x).GetBondTypeAsDouble() == 2.0:
            decorations.append(("oxo", x, None))
            continue
        if mol.GetBondBetweenAtoms(carbon, x).GetBondTypeAsDouble() != 1.0:
            return None
        if z == 8 and atom.GetDegree() == 1:
            decorations.append(("OH", x, None))
        elif z == 8 and atom.GetDegree() == 2 and _acyl_carbon(mol, graph, x, carbon, piece):
            decorations.append(("acyl", x, next(n for n in graph[x] if n != carbon)))
        elif z == 8 and atom.GetDegree() == 2:
            other = next(n for n in graph[x] if n != carbon)
            if mol.GetAtomWithIdx(other).GetAtomicNum() != 6 or other in piece:
                return None
            if any(mol.GetAtomWithIdx(a).GetAtomicNum() not in (6, 8, 9, 17, 35, 53) for a in subtree_atoms(graph, other, x)):
                return None
            decorations.append(("ether", x, other))
        elif z == 7:
            decorations.append(("N", x, None))
        elif z in _CHAIN_HALOGENS and atom.GetDegree() == 1:
            decorations.append(("X", x, None))
        else:
            return None
    covered = set(piece)
    for (kind, x, root), carbon in zip(decorations, chain):
        if kind == "COOH":
            covered |= {n for n in graph[carbon] if n not in piece}
        elif x is not None:
            covered |= subtree_atoms(graph, x, carbon)
    return (chain, decorations) if covered == set(range(mol.GetNumAtoms())) else None


def _acyl_carbon(mol, graph, oxygen, carbon, piece):
    other = next((n for n in graph[oxygen] if n != carbon), None)
    if other is None or other in piece or mol.GetAtomWithIdx(other).GetAtomicNum() != 6:
        return False
    return any(
        mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(other, n).GetBondTypeAsDouble() == 2.0
        for n in graph[other]
    )


def subtree_atoms(graph, root, blocked):
    from ._cited_group import subtree

    return subtree(graph, root, blocked)


def _restored_chain(mol, chain, decorations):
    """The plain polyhydroxy chain (every position CH-OH, both ends CH2OH or the carboxy end kept) with its atoms."""
    from ._cited_group import subtree

    graph = adjacency(mol)
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    removed = set()
    for carbon, (kind, x, root) in zip(chain, decorations):
        if kind == "H":
            oxygen = editable.AddAtom(Chem.Atom(8))
            editable.GetAtomWithIdx(oxygen).SetIntProp("_orig", -1)
            editable.AddBond(carbon, oxygen, Chem.BondType.SINGLE)
        elif kind in ("N", "X"):
            editable.GetAtomWithIdx(x).SetAtomicNum(8)
            for n in graph[x]:
                if n != carbon:
                    removed |= subtree(graph, n, x)
        elif kind in ("ether", "acyl"):
            removed |= subtree(graph, root, x)
        elif kind == "oxo":
            editable.GetBondBetweenAtoms(carbon, x).SetBondType(Chem.BondType.SINGLE)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    model = editable.GetMol()
    Chem.SanitizeMol(model)
    return model, {a.GetIntProp("_orig"): a.GetIdx() for a in model.GetAtoms() if a.GetIntProp("_orig") >= 0}


def substituted_chain_name(mol, bridges=()):
    """P-102.5.6.5, P-102.5.6.6.2: an alditol or aldonic acid with deoxy, amino, halogen or ether positions: the chain
    is restored to the plain alditol, named from its aldose, and the substituents are cited as prefixes."""
    from ._cited_group import cited_group
    from ._sugar_substituted import _HALOGENS, _plain_group, _segment

    found = _decorated_chain(mol)
    if found is None:
        return None
    chain, decorations = found
    if all(kind in ("OH", "COOH") for kind, _, _ in decorations) and not bridges:
        return None
    if sum(kind in ("OH", "ether", "acyl") for kind, _, _ in decorations) < 3:
        return None
    try:
        model, atom_of = _restored_chain(mol, chain, decorations)
    except (ValueError, RuntimeError, Chem.rdchem.MolSanitizeException):
        return None
    modeled = _chain(model)
    if modeled is None:
        return None
    model_chain = [atom_of[c] for c in chain]
    if model_chain != modeled[0] and model_chain != list(reversed(modeled[0])):
        return None
    kinds = modeled[1] if model_chain == modeled[0] else list(reversed(modeled[1]))
    kind = _class(kinds)
    if kind not in ("alditol", "aldonic"):
        return None
    graph = adjacency(mol)
    candidates = []
    for reverse in ([False, True] if kind == "alditol" else [kinds[0][0] != "COOH"]):
        order = list(reversed(range(len(chain)))) if reverse else list(range(len(chain)))
        try:
            aldose = _as_aldose(model, model_chain, kinds, reverse)
        except Exception:
            aldose = None
        if aldose is None or not aldose.endswith("ose"):
            continue
        entries, esters, bridged, amino = {}, {}, [], []
        failed = False
        for position, index in enumerate(order, start=1):
            dkind, x, root = decorations[index]
            if dkind in ("OH", "COOH"):
                continue
            if dkind == "acyl" and kind == "alditol":
                from ._inositol_derivative import ester_anion

                anion = ester_anion(mol, x, root)
                if anion is None:
                    failed = True
                    break
                esters.setdefault(anion, []).append(position)
                continue
            if dkind == "acyl":
                name, compound = cited_group(mol, graph, root, x)
                entries.setdefault((name, name, "O", compound), []).append(position)
                continue
            if dkind == "H":
                entries.setdefault(("deoxy", "deoxy", "", False), []).append(position)
            elif dkind == "N":
                name, compound = cited_group(mol, graph, x, chain[index])
                entries.setdefault((name, name, "", compound), []).append(position)
                entries.setdefault(("deoxy", "deoxy", "", False), []).append(position)
                amino.append((position, name, compound, (mol, graph, x, chain[index])))
            elif dkind == "X":
                text = _HALOGENS[mol.GetAtomWithIdx(x).GetAtomicNum()]
                entries.setdefault((text, text, "", False), []).append(position)
                entries.setdefault(("deoxy", "deoxy", "", False), []).append(position)
            elif dkind == "ether":
                if not _plain_group(mol, graph, root, x):
                    failed = True
                    break
                name, compound = cited_group(mol, graph, root, x)
                entries.setdefault((name, name, "O", compound), []).append(position)
        if failed:
            continue
        for first, second in bridges:
            if first not in chain or second not in chain:
                failed = True
                break
            ends = sorted((order.index(chain.index(first)) + 1, order.index(chain.index(second)) + 1))
            bridged.append(ends)
        if failed:
            continue
        candidates.append((aldose, entries, esters, bridged, amino))
    if not candidates:
        return None
    from ._common import alpha_sort_key

    def rank(item):
        aldose, entries, esters, bridged, _ = item
        locants = sorted([p for ps in [*entries.values(), *esters.values()] for p in ps] + [p for ends in bridged for p in ends])
        first = min(entries, key=lambda key: alpha_sort_key(key[0])) if entries else None
        return (aldose[2:], aldose[0] != "D", -len(locants), locants, min(entries[first]) if first else 0)

    if kind == "alditol":
        for aldose, entries, esters, bridged, amino in candidates:
            if len(amino) == 1 and amino[0][0] == 2 and not esters and not bridged and aldose in _AMINO_ALDITOLS:
                from ._sugar_substituted import glycosamine_name

                deoxy = ("deoxy", "deoxy", "", False)
                retained = glycosamine_name(aldose, entries.get(deoxy, []), amino[0][3])
                if retained is not None and len(entries) == 2:
                    prefix = "".join(text + "-" for _, text in sorted(retained[1], key=lambda item: (alpha_sort_key(item[0]), item[1])))
                    return prefix + retained[0][: -len("e")] + "itol"
    if kind == "alditol":
        for aldose, entries, esters, bridged, _ in candidates:
            if not esters and not bridged and list(entries) == [("deoxy", "deoxy", "", False)] and entries[("deoxy", "deoxy", "", False)] == [len(chain)] and aldose in _RETAINED_ALDITOLS:
                return _RETAINED_ALDITOLS[aldose]
    aldose, entries, esters, bridged, _ = min(candidates, key=rank)
    ending = _suffixes(kind)
    stem = aldose[: -len("ose")]
    cited = [(key[0], _segment(sorted(locants), key[0], key[2], key[3])) for key, locants in entries.items()]
    if bridged:
        from ._sugar_substituted import anhydro_segment

        cited.append(("anhydro", anhydro_segment(bridged)))
    segments = [text for _, text in sorted(cited, key=lambda item: (alpha_sort_key(item[0]), item[1]))]
    prefix = "-".join(segments) + "-" if segments else ""
    name = prefix + stem + ending
    if esters:
        from ._inositol_derivative import ester_words

        words = ester_words([(a, sorted(l)) for a, l in esters.items()])
        if len(esters) == 1 and sum(map(len, esters.values())) == len(chain) and not entries:
            words = words.split("-", 1)[1]
        name = f"{name} {words}"
    return name


_AMINO_ALDITOLS = {"D-glucose", "D-galactose"}


def has_substituted_chain_shape(mol) -> bool:
    try:
        return substituted_chain_name(mol) is not None
    except (UnsupportedStructure, ValueError, RuntimeError):
        return False


def name_substituted_chain(mol) -> str:
    return substituted_chain_name(mol)


_RING_AMIDE = Chem.MolFromSmarts("[CX3;R;H0](=O)[NX3;R;H1][CX4;R]")


def sugar_lactam_name(mol):
    """'5-amino-5-deoxy-D-galactono-1,5-lactam' for the internal amide of an amino aldonic acid (P-102.5.6.6.2.2): the
    ring is opened at the amide nitrogen, the amino acid named, and the locants of the carboxy carbon and of the
    carbon that carries the nitrogen cited."""
    matches = mol.GetSubstructMatches(_RING_AMIDE)
    if len(matches) != 1 or mol.GetRingInfo().NumRings() != 1:
        return None
    carbonyl, _, nitrogen, closing = matches[0]
    editable = Chem.RWMol(mol)
    editable.RemoveBond(carbonyl, nitrogen)
    hydroxy = editable.AddAtom(Chem.Atom(8))
    editable.AddBond(carbonyl, hydroxy, Chem.BondType.SINGLE)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    try:
        chain_mol = editable.GetMol()
        Chem.SanitizeMol(chain_mol)
        acid = substituted_chain_name(chain_mol)
    except (ValueError, RuntimeError, Chem.rdchem.MolSanitizeException):
        return None
    if acid is None or not acid.endswith("onic acid"):
        return None
    position, seen, frontier = {carbonyl: 1}, {carbonyl}, [carbonyl]
    while frontier:
        atom = frontier.pop(0)
        for n in chain_mol.GetAtomWithIdx(atom).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIdx() not in seen:
                seen.add(n.GetIdx())
                position[n.GetIdx()] = position[atom] + 1
                frontier.append(n.GetIdx())
    return f"{acid[: -len('ic acid')]}o-1,{position[closing]}-lactam"


def has_sugar_lactam_shape(mol) -> bool:
    try:
        return sugar_lactam_name(mol) is not None
    except (UnsupportedStructure, ValueError, RuntimeError):
        return False


def name_sugar_lactam(mol) -> str:
    return sugar_lactam_name(mol)


_RING_ACID = Chem.MolFromSmarts("[CX4;R]([OX2;R])-[CX3;H0;!R](=O)[OX2,NX3]")


def _nitrogen_prefix(substituents):
    names = sorted(name for name, _ in substituents)
    if not names:
        return ""
    if len(names) == 1:
        return f"N-{names[0]}"
    if len(set(names)) == 1:
        return f"N,N-di{names[0]}" if not any(ch in names[0] for ch in "-()[], ") else f"N,N-di({names[0]})"
    return None


def sugar_ring_acid_name(mol):
    """P-102.5.6.6.3, P-102.5.6.6.4: ring-closed uronic acids and ketoaldonic acids and their derivatives. The uronic
    derivative is turned into the acid, the carboxy group of a ketoaldonic acid into the CH2OH of the ketose; the sugar
    is named and its ending becomes 'uronic acid', 'onic acid' or the ending of the ester or amide, a glycoside being
    isolated in parentheses when it is esterified or amidated."""
    from ._cited_group import cited_group, subtree
    from .core import smiles_to_iupac

    matches = mol.GetSubstructMatches(_RING_ACID)
    if len(matches) != 1 or mol.GetRingInfo().NumRings() != 1:
        return None
    ring_carbon, ring_oxygen, carboxy, oxo, hetero = matches[0]
    graph = adjacency(mol)
    ketose = any(
        mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and n != ring_oxygen and mol.GetAtomWithIdx(n).GetDegree() <= 2
        for n in graph[ring_carbon]
        if n != carboxy
    )
    het = mol.GetAtomWithIdx(hetero)
    removal, alcohol, substituents = set(), None, []
    if het.GetAtomicNum() == 8 and het.GetDegree() == 1:
        derivative = "acid"
    elif het.GetAtomicNum() == 8:
        derivative = "ester"
        root = next(n for n in graph[hetero] if n != carboxy)
        alcohol = cited_group(mol, graph, root, hetero)[0]
        removal |= subtree(graph, root, hetero)
    else:
        derivative = "amide"
        for n in graph[hetero]:
            if n != carboxy:
                substituents.append(cited_group(mol, graph, n, hetero))
                removal |= subtree(graph, n, hetero)
    if derivative == "acid" and not ketose:
        return None
    editable = Chem.RWMol(mol)
    if ketose:
        removal |= {oxo, hetero}
        hydroxy = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(carboxy, hydroxy, Chem.BondType.SINGLE)
    elif derivative == "amide":
        editable.GetAtomWithIdx(hetero).SetAtomicNum(8)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    for index in sorted(removal, reverse=True):
        editable.RemoveAtom(index)
    try:
        model = editable.GetMol()
        Chem.SanitizeMol(model)
        name = smiles_to_iupac(Chem.MolToSmiles(model))
    except (UnsupportedStructure, ValueError, RuntimeError, Chem.rdchem.MolSanitizeException):
        return None
    nitrogen = _nitrogen_prefix(substituents)
    if nitrogen is None:
        return None
    glycoside = name.endswith("oside") or "osiduronic acid" in name
    if ketose:
        if not name.endswith(("ose", "oside")):
            return None
        stem = name[:-1]
        if derivative == "acid":
            return f"{stem}onic acid"
        if derivative == "ester":
            return f"{alcohol} ({stem})onate" if glycoside else f"{alcohol} {stem}onate"
        return f"{nitrogen}({stem})onamide" if glycoside else f"{nitrogen + '-' if nitrogen else ''}{stem}onamide"
    if not name.endswith("uronic acid"):
        return None
    if derivative == "acid":
        return name
    stem = name[: -len("ic acid")]
    base = name[: -len("uronic acid")]
    if derivative == "ester":
        return f"{alcohol} ({base})uronate" if glycoside else f"{alcohol} {stem}ate"
    if glycoside and nitrogen:
        return f"{nitrogen}({base})uronamide"
    return f"{nitrogen + '-' if nitrogen else ''}{stem}amide"


def has_sugar_ring_acid_shape(mol) -> bool:
    try:
        return sugar_ring_acid_name(mol) is not None
    except (UnsupportedStructure, ValueError, RuntimeError):
        return False


def name_sugar_ring_acid(mol) -> str:
    return sugar_ring_acid_name(mol)


def ketoaldonic_chain_name(mol):
    """P-102.5.6.6.3.1: an open-chain ketoaldonic acid, 'D-arabino-hex-5-ulosonic acid': the keto group is restored to a
    hydroxy group for the configurational prefixes of the centres that remain, numbering starts at the carboxy group."""
    from ._cited_group import cited_group
    from ._common import alpha_sort_key
    from ._sugar_substituted import _HALOGENS, _Skeleton, _labels, _parent, _plain_group, _segment

    found = _decorated_chain(mol, allow_oxo=True)
    if found is None:
        return None
    chain, decorations = found
    if sum(kind == "oxo" for kind, _, _ in decorations) != 1 or sum(kind == "COOH" for kind, _, _ in decorations) != 1:
        return None
    if decorations[-1][0] == "COOH":
        chain, decorations = list(reversed(chain)), list(reversed(decorations))
    if decorations[0][0] != "COOH" or decorations[-1][0] == "oxo":
        return None
    if sum(kind in ("OH", "ether", "acyl") for kind, _, _ in decorations) < 2:
        return None
    graph = adjacency(mol)
    try:
        model, atom_of = _restored_chain(mol, chain, decorations)
        skeleton = _Skeleton(chain, False)
        centres, _ = _labels(model, skeleton, atom_of)
        parent = _parent(skeleton, centres, None, False)
    except (UnsupportedStructure, ValueError, RuntimeError, KeyError, Chem.rdchem.MolSanitizeException):
        return None
    if parent is None:
        return None
    core = parent[0]
    if not core.endswith("ose"):
        return None
    position = next(i for i, (kind, _, _) in enumerate(decorations, start=1) if kind == "oxo")
    entries = {}
    for place, (kind, x, root) in enumerate(decorations, start=1):
        if kind in ("OH", "COOH", "oxo"):
            continue
        carbon = chain[place - 1]
        if kind == "H":
            entries.setdefault(("deoxy", "deoxy", "", False), []).append(place)
        elif kind == "N":
            name, compound = cited_group(mol, graph, x, carbon)
            entries.setdefault((name, name, "", compound), []).append(place)
            entries.setdefault(("deoxy", "deoxy", "", False), []).append(place)
        elif kind == "X":
            text = _HALOGENS[mol.GetAtomWithIdx(x).GetAtomicNum()]
            entries.setdefault((text, text, "", False), []).append(place)
            entries.setdefault(("deoxy", "deoxy", "", False), []).append(place)
        elif kind == "ether":
            if not _plain_group(mol, graph, root, x):
                return None
            name, compound = cited_group(mol, graph, root, x)
            entries.setdefault((name, name, "O", compound), []).append(place)
        elif kind == "acyl":
            name, compound = cited_group(mol, graph, root, x)
            entries.setdefault((name, name, "O", compound), []).append(place)
    cited = [(key[0], _segment(sorted(locants), key[0], key[2], key[3])) for key, locants in entries.items()]
    segments = [text for _, text in sorted(cited, key=lambda item: (alpha_sort_key(item[0]), item[1]))]
    prefix = "-".join(segments) + "-" if segments else ""
    return f"{prefix}{core[: -len('ose')]}-{position}-ulosonic acid"


def has_ketoaldonic_chain_shape(mol) -> bool:
    try:
        return ketoaldonic_chain_name(mol) is not None
    except (UnsupportedStructure, ValueError, RuntimeError):
        return False


def name_ketoaldonic_chain(mol) -> str:
    return ketoaldonic_chain_name(mol)
