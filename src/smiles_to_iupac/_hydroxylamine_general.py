"""O-substituted hydroxylamines and their chalcogen analogues, single or multiplied (P-68.3.1.1.1.2, P-68.3.1.1.1.6):
H2N-O-R is 'O-(chloromethyl)hydroxylamine', H2N-S-R is 'S-methyl(thiohydroxylamine)', and two or more such groups on one
carbon skeleton are 'O,O'-(ethane-1,2-diyl)bis(hydroxylamine)'. The group R is named by the substituent machinery."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_PARENT = {8: ("O", "hydroxylamine"), 16: ("S", "thiohydroxylamine"), 34: ("Se", "selenohydroxylamine"), 52: ("Te", "tellurohydroxylamine")}


def _groups(mol):
    """[(nitrogen, chalcogen, carbon)] of every H2N-X-C group, or None when the molecule is anything else."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    found = []
    for nitrogen in (a for a in mol.GetAtoms() if a.GetAtomicNum() == 7):
        if nitrogen.GetDegree() != 1 or nitrogen.GetTotalNumHs() != 2 or nitrogen.GetFormalCharge():
            return None
        (chalcogen,) = nitrogen.GetNeighbors()
        if chalcogen.GetAtomicNum() not in _PARENT or chalcogen.GetDegree() != 2 or chalcogen.GetFormalCharge() or chalcogen.IsInRing():
            return None
        root = next(n for n in chalcogen.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx())
        if root.GetAtomicNum() not in (6, 16) or (root.GetAtomicNum() == 16 and chalcogen.GetAtomicNum() != 8):
            return None
        found.append((nitrogen, chalcogen, root))
    if not found or len({c.GetAtomicNum() for _, c, _ in found}) != 1:
        return None
    members = {a.GetIdx() for n, c, _ in found for a in (n, c)}
    for _, _, root in found:
        if root.GetAtomicNum() == 16:
            oxo = [n for n in root.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
            members.update([root.GetIdx(), *(n.GetIdx() for n in oxo)])
    if any(
        a.GetIdx() not in members and (a.GetAtomicNum() not in (6, 9, 17, 35, 53) or a.GetFormalCharge() or a.GetIsotope())
        for a in mol.GetAtoms()
    ):
        return None
    return found


def _bare_chalcogen_analogue(mol):
    """'thiohydroxylamine' for H2N-SH and its selenium and tellurium analogues (the oxygen parent is hydroxylamine)."""
    if mol.GetNumAtoms() != 2 or Chem.GetMolFrags(mol).__len__() != 1:
        return None
    nitrogen = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    chalcogen = [a for a in mol.GetAtoms() if a.GetAtomicNum() in (16, 34, 52)]
    if len(nitrogen) == 1 and len(chalcogen) == 1 and not any(a.GetFormalCharge() for a in mol.GetAtoms()):
        return _PARENT[chalcogen[0].GetAtomicNum()][1]
    return None


def has_o_substituted_hydroxylamine_shape(mol) -> bool:
    return _groups(mol) is not None or _bare_chalcogen_analogue(mol) is not None


def name_o_substituted_hydroxylamine(mol) -> str:
    found = _groups(mol)
    if found is None:
        bare = _bare_chalcogen_analogue(mol)
        if bare is None:
            raise UnsupportedStructure("not an O-substituted hydroxylamine")
        return bare
    symbol, parent = _PARENT[found[0][1].GetAtomicNum()]
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    if len(found) == 1:
        nitrogen, chalcogen, root = found[0]
        name, compound = name_branch(graph, root.GetIdx(), chalcogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
        prefix = format_substituent_prefixes({name: {"locants": [symbol], "compound": compound}})
        return prefix + (parent if symbol == "O" else f"({parent})")
    return _multiplied_name(mol, graph, found, symbol, parent)


def _multiplied_name(mol, graph, found, symbol, parent):
    from ._diester_ring_diyl import evaluate_skeleton, select_skeleton

    matches = [(None, None, chalcogen, root) for _, chalcogen, root in found]
    selection = select_skeleton(mol, graph, matches)
    if selection is None:
        raise UnsupportedStructure("no skeleton carries the hydroxylamine groups")
    kind, body, pool = selection
    options = [(kind, [path], set(path)) for path in body] if kind == "chain" else [(kind, body, pool)]
    attach = [root.GetIdx() for _, _, root in found]
    blocked = {chalcogen.GetIdx() for _, chalcogen, _ in found}
    best = None
    for option_kind, option_body, option_pool in options:
        result = evaluate_skeleton(mol, graph, option_kind, option_body, option_pool, attach, blocked, "yl")
        if result is not None and (best is None or result[0] < best[0]):
            best = result
    if best is None:
        raise UnsupportedStructure("no admissible numbering places every hydroxylamine group on this skeleton")
    locants = ",".join(symbol + "'" * k for k in range(len(found)))
    return f"{locants}-({best[1]}){multiplying_prefix(len(found), compound=True)}({parent})"
