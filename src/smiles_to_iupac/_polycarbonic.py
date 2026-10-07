"""Di-, tri-, tetra- and polycarbonic acids and their functional-replacement analogues (P-65.2.3): HO-[CO-O]n-H numbered
along the chain, carbons at odd positions; replacement is cited by prefixes with locants, omitted when every oxygen is
replaced by one prefix (P-65.2.3.1.2.1); five or more units are named by skeletal replacement.
"""

from rdkit import Chem

from ._acid_lexicon import _INFIX
from ._carbonic_family import _HALIDE_INFIX, _SYMBOL, PSEUDO_INFIX, _chain, pseudohalide_at
from ._common import UnsupportedStructure, adjacency
from ._numerals import multiplying_prefix

_PREFIX_WORDS = {"fluorid": "fluoro", "chlorid": "chloro", "bromid": "bromo", "iodid": "iodo", "azid": "azido", "isocyanid": "isocyano", "isocyanatid": "isocyanato", "isothiocyanatid": "isothiocyanato", "isoselenocyanatid": "isoselenocyanato", "isotellurocyanatid": "isotellurocyanato"}
_UNITS = {2: "di", 3: "tri", 4: "tetra"}
_BRIDGE_WORD = {"S": "thio", "Se": "seleno", "Te": "telluro", "NH": "imido"}


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _carbonic_center(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetFormalCharge() or atom.GetTotalNumHs() or atom.GetDegree() != 3:
        return False
    return sum(1 for n in atom.GetNeighbors() if _bond(mol, idx, n.GetIdx()) == 2.0) == 1


def _oxo_word(mol, center):
    for n in mol.GetAtomWithIdx(center).GetNeighbors():
        if _bond(mol, center, n.GetIdx()) == 2.0:
            z = n.GetAtomicNum()
            if z in _SYMBOL and n.GetDegree() == 1:
                return _SYMBOL[z], n.GetIdx()
            if z == 7 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1:
                return "NH", n.GetIdx()
    return None, None


def name_polycarbonic(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("several fragments")
    centers = [a.GetIdx() for a in mol.GetAtoms() if _carbonic_center(mol, a.GetIdx())]
    if len(centers) < 2:
        raise UnsupportedStructure("fewer than two carbonic units")
    graph = adjacency(mol)
    bridges = {}
    for c in centers:
        for n in graph[c]:
            atom = mol.GetAtomWithIdx(n)
            if atom.GetAtomicNum() in _SYMBOL and atom.GetDegree() == 2 and not atom.GetFormalCharge() and not atom.GetTotalNumHs():
                far = [m for m in graph[n] if m != c]
                if far and far[0] in centers and _bond(mol, c, n) == 1.0 and _bond(mol, n, far[0]) == 1.0:
                    bridges.setdefault(frozenset((c, far[0])), n)
            elif atom.GetAtomicNum() == 7 and atom.GetDegree() == 2 and atom.GetTotalNumHs() == 1 and not atom.GetFormalCharge():
                far = [m for m in graph[n] if m != c]
                if far and far[0] in centers:
                    bridges.setdefault(frozenset((c, far[0])), n)
    if len(bridges) != len(centers) - 1:
        raise UnsupportedStructure("the carbonic units are not joined in a single chain")
    degree = {c: sum(1 for pair in bridges if c in pair) for c in centers}
    ends = [c for c in centers if degree[c] == 1]
    if len(ends) != 2 or any(d > 2 for d in degree.values()):
        raise UnsupportedStructure("a branched polycarbonic acid is not named by P-65.2.3")
    n = len(centers)
    best = None
    for start in ends:
        order = [start]
        while len(order) < n:
            order.append(next(m for pair in bridges if order[-1] in pair for m in pair if m != order[-1] and m not in order))
        replaced = _replacements(mol, graph, order, bridges)
        key = tuple(sorted(loc for loc, _ in replaced))
        if best is None or key < best[0]:
            best = (key, order, replaced)
    _, order, replaced = best
    ends_ligands = [
        w for _, w in replaced if w in _TERMINAL_ONLY
    ]
    if len(ends_ligands) == 2:
        raise UnsupportedStructure("two non-hydroxy end groups make an acyl halide of the polycarbonic acid")
    return _compose(mol, n, replaced, _letters(mol, graph, order, bridges))


def _replacements(mol, graph, order, bridges):
    """[(locant, prefix word)] and letters for the replaced atoms along the numbered chain."""
    n = len(order)
    replaced, letters = [], []
    for i, center in enumerate(order):
        position = 2 * i + 1
        oxo, _ = _oxo_word(mol, center)
        if oxo is None:
            raise UnsupportedStructure("unsupported =X atom on a polycarbonic acid")
        if oxo != "O":
            replaced.append((position, _BRIDGE_WORD.get(oxo) or {"NNH2": "hydrazono"}.get(oxo, _INFIX.get(oxo))))
        if i in (0, n - 1):
            bridge_atoms = set(bridges.values())
            ligand = [m for m in graph[center] if m not in bridge_atoms and _bond(mol, center, m) == 1.0]
            if len(ligand) != 1:
                raise UnsupportedStructure("a terminal group of the polycarbonic acid is not recognized")
            _terminal(mol, graph, center, ligand[0], position, replaced, letters)
    for i in range(n - 1):
        bridge = bridges[frozenset((order[i], order[i + 1]))]
        atom = mol.GetAtomWithIdx(bridge)
        symbol = "NH" if atom.GetAtomicNum() == 7 else _SYMBOL[atom.GetAtomicNum()]
        if symbol != "O":
            replaced.append((2 * i + 2, _BRIDGE_WORD[symbol]))
    return replaced


def _terminal(mol, graph, center, atom_idx, position, replaced, letters):
    atom = mol.GetAtomWithIdx(atom_idx)
    chain = _chain(mol, center, atom)
    if chain is not None:
        symbols = chain[0]
        if len(symbols) == 2:
            if symbols != ["O", "O"]:
                raise UnsupportedStructure("a mixed peroxy chalcogen group is not handled on a polycarbonic acid")
            replaced.append((position, "peroxy"))
        elif symbols[0] != "O":
            replaced.append((position, _INFIX[symbols[0]]))
        return
    if atom.GetAtomicNum() in _HALIDE_INFIX and atom.GetDegree() == 1:
        replaced.append((position, _PREFIX_WORDS[_HALIDE_INFIX[atom.GetAtomicNum()]]))
        return
    pseudo = pseudohalide_at(mol, atom_idx, center)
    if pseudo is not None and pseudo != "CN":
        replaced.append((position, _PREFIX_WORDS[PSEUDO_INFIX[pseudo]]))
        return
    raise UnsupportedStructure("this terminal group of a polycarbonic acid is not supported")


_TERMINAL_ONLY = {"fluoro", "chloro", "bromo", "iodo", "azido", "isocyano", "isocyanato", "isothiocyanato", "isoselenocyanato", "isotellurocyanato"}


def _letters(mol, graph, order, bridges):
    """Italic letters with locants of the hydrogen-bearing atoms of the two end groups when a replaced
    chalcogen could sit on either the =X or the -Y position of that carbon (P-65.2.3.1.2.2)."""
    n = len(order)
    found = []
    bridge_atoms = set(bridges.values())
    for i in (0, n - 1):
        center = order[i]
        x, _ = _oxo_word(mol, center)
        ligand = [m for m in graph[center] if m not in bridge_atoms and _bond(mol, center, m) == 1.0]
        chain = _chain(mol, center, mol.GetAtomWithIdx(ligand[0])) if ligand else None
        if chain is None or len(chain[0]) != 1:
            continue
        y = chain[0][0]
        if (x != "O") != (y != "O") and x in ("O", "S", "Se", "Te") :
            found.append((f"{y}{2 * i + 1}", y != "O"))
    if not found:
        return ""
    marked = [text for text, non_oxygen in found if non_oxygen] or [text for text, _ in found]
    return ",".join(sorted(marked))


def _compose(mol, n, replaced, letters=""):
    if n >= 5:
        if replaced:
            raise UnsupportedStructure("functional replacement on a polycarbonic acid with five or more units is not named")
        interior = ",".join(str(p) for p in range(3, 2 * n - 2, 2))
        oxa = ",".join(str(p) for p in range(2, 2 * n - 1, 2))
        length = {9: "nonane", 11: "undecane", 13: "tridecane", 15: "pentadecane", 17: "heptadecane", 19: "nonadecane"}.get(2 * n - 1)
        if length is None:
            raise UnsupportedStructure("this polycarbonic acid length has no stem")
        return (
            f"{interior}-{multiplying_prefix(n - 2)}oxo-{oxa}-{multiplying_prefix(n - 1)}oxa{length[:-1]}edioic acid"
        )
    core = _UNITS[n] + "carbonic acid"
    if not replaced:
        return core
    by_word = {}
    for position, word in replaced:
        by_word.setdefault(word, []).append(position)
    total_oxygens = 2 * n + 1
    parts = []
    for word in sorted(by_word):
        positions = sorted(by_word[word])
        cited = "" if len(by_word) == 1 and word in ("thio", "seleno", "telluro") and len(positions) == total_oxygens else ",".join(map(str, positions)) + "-"
        multiplier = multiplying_prefix(len(positions)) if len(positions) > 1 else ""
        parts.append(f"{cited}{multiplier}{word}")
    only_terminal = all(word in _TERMINAL_ONLY for word in by_word)
    if only_terminal:
        parts = [
            (multiplying_prefix(len(by_word[w])) if len(by_word[w]) > 1 else "") + w for w in sorted(by_word)
        ]
    return "-".join(parts) + core[:-len(" acid")] + (f" {letters}-acid" if letters else " acid")
