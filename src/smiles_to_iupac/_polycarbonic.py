"""Di-, tri-, tetra- and polycarbonic acids and their functional-replacement analogues (P-65.2.3): HO-[CO-O]n-H numbered
along the chain, carbons at odd positions; replacement is cited by prefixes with locants, omitted when every oxygen is
replaced by one prefix (P-65.2.3.1.2.1); five or more units are named by skeletal replacement.
"""

import re

from rdkit import Chem

from ._acid_lexicon import _INFIX
from ._carbonic_family import _HALIDE_INFIX, _SYMBOL, PSEUDO_INFIX, _chain, pseudohalide_at
from ._alkoxy import alkoxy_prefix
from ._carbonic_amide import _arm_atoms
from ._common import UnsupportedStructure, adjacency, halogen_substituents, is_halide_ligand
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._polycarbonic_nuclides import chain_nuclide_key, chain_nuclides_text, insert_chain_nuclides
from ._substituents import ISOTOPE_LABELS, format_substituent_prefixes, name_branch

_PREFIX_WORDS = {"fluorid": "fluoro", "chlorid": "chloro", "bromid": "bromo", "iodid": "iodo", "azid": "azido", "isocyanid": "isocyano", "isocyanatid": "isocyanato", "isothiocyanatid": "isothiocyanato", "isoselenocyanatid": "isoselenocyanato", "isotellurocyanatid": "isotellurocyanato"}
_UNITS = {2: "di", 3: "tri", 4: "tetra"}
_BRIDGE_WORD = {"S": "thio", "Se": "seleno", "Te": "telluro", "NH": "imido"}


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _carbonic_center(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetFormalCharge() or atom.GetTotalNumHs() or atom.GetDegree() != 3:
        return False
    if any(n.GetAtomicNum() == 6 and pseudohalide_at(mol, n.GetIdx(), idx) != "CN" for n in atom.GetNeighbors()):
        return False
    return sum(1 for n in atom.GetNeighbors() if _bond(mol, idx, n.GetIdx()) == 2.0) == 1


def _oxo_word(mol, center, substituted=False):
    for n in mol.GetAtomWithIdx(center).GetNeighbors():
        if _bond(mol, center, n.GetIdx()) == 2.0:
            z = n.GetAtomicNum()
            if z in _SYMBOL and n.GetDegree() == 1:
                return _SYMBOL[z], n.GetIdx()
            if z == 7 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1:
                return "NH", n.GetIdx()
            if z == 7 and n.GetDegree() == 2 and not n.GetTotalNumHs() and not n.GetFormalCharge() and not n.IsInRing():
                far = next(m for m in n.GetNeighbors() if m.GetIdx() != center)
                if far.GetAtomicNum() == 7:
                    if far.GetDegree() == 1 and far.GetTotalNumHs() == 2 and not far.GetFormalCharge():
                        return "NNH2", n.GetIdx()
                    return None, None
                if substituted:
                    return "NH", n.GetIdx()
    return None, None


def _dianion_name(mol):
    """'1-thiodicarbonate': the name of the fully deprotonated acid, whose tautomer letters are no longer needed."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 2 or any(a.GetFormalCharge() != -1 or a.GetDegree() != 1 or a.GetAtomicNum() not in _SYMBOL for a in charged):
        raise UnsupportedStructure("only the fully deprotonated polycarbonic acids are named as anions here")
    neutral = Chem.RWMol(mol)
    for atom in charged:
        target = neutral.GetAtomWithIdx(atom.GetIdx())
        target.SetFormalCharge(0)
        target.SetNoImplicit(False)
    neutral = neutral.GetMol()
    Chem.SanitizeMol(neutral)
    name = name_polycarbonic(neutral)
    name = re.sub(r" (?:[A-Za-z]+\d+)(?:,[A-Za-z]+\d+)*-acid$", " acid", name)
    return name[: -len("ic acid")] + "ate"


def _chain_setup(mol, substituted=False):
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
            if atom.GetAtomicNum() == 8 and atom.GetDegree() == 2 and not atom.GetFormalCharge() and not atom.GetTotalNumHs():
                partner = [m for m in graph[n] if m != c]
                if partner and mol.GetAtomWithIdx(partner[0]).GetAtomicNum() == 8 and mol.GetAtomWithIdx(partner[0]).GetDegree() == 2:
                    far = [m for m in graph[partner[0]] if m != n]
                    if far and far[0] in centers:
                        bridges.setdefault(frozenset((c, far[0])), (n, partner[0]))
            elif atom.GetAtomicNum() == 7 and not atom.GetFormalCharge() and not atom.IsInRing() and (
                (atom.GetDegree() == 2 and atom.GetTotalNumHs() == 1) or (substituted and atom.GetDegree() == 3 and not atom.GetTotalNumHs())
            ):
                far = [m for m in graph[n] if m != c and m in centers]
                if len(far) == 1:
                    bridges.setdefault(frozenset((c, far[0])), n)
    if len(bridges) != len(centers) - 1:
        raise UnsupportedStructure("the carbonic units are not joined in a single chain")
    degree = {c: sum(1 for pair in bridges if c in pair) for c in centers}
    ends = [c for c in centers if degree[c] == 1]
    if len(ends) != 2 or any(d > 2 for d in degree.values()):
        raise UnsupportedStructure("a branched polycarbonic acid is not named by P-65.2.3")
    return graph, bridges, ends, len(centers)


def _walk(bridges, start, n):
    order = [start]
    while len(order) < n:
        order.append(next(m for pair in bridges if order[-1] in pair for m in pair if m != order[-1] and m not in order))
    return order


def end_numberings(mol):
    """P-65.2.3.1: the chain orders of a di- to tetracarbonic acid that give the replacement prefixes the lowest
    locants, as [(order, locant of each carbonic carbon)]; the carbons sit at the odd positions."""
    graph, bridges, ends, n = _chain_setup(mol)
    if n > 4:
        raise UnsupportedStructure("five or more carbonic units are named by skeletal replacement")
    keyed = []
    for start in ends:
        order = _walk(bridges, start, n)
        replaced = _replacements(mol, graph, order, bridges)
        keyed.append((tuple(sorted(loc for loc, _ in replaced)), order))
    lowest = min(key for key, _ in keyed)
    return [(order, {c: 2 * i + 1 for i, c in enumerate(order)}) for key, order in keyed if key == lowest]


def name_polycarbonic(mol, numbered=None):
    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("several fragments")
    if any(a.GetFormalCharge() for a in mol.GetAtoms()):
        return _dianion_name(mol)
    graph, bridges, ends, n = _chain_setup(mol, substituted=True)
    found = []
    for start in ends:
        order = _walk(bridges, start, n)
        try:
            replaced = _replacements(mol, graph, order, bridges, substituted=True)
        except UnsupportedStructure:
            acyl = _acyl_end_name(mol, graph, order, bridges)
            if acyl is None:
                raise
            return acyl
        entries = _nitrogen_substituents(mol, graph, order, bridges)
        key = (
            tuple(sorted(loc for loc, _ in replaced)),
            tuple(sorted(_locant_number(loc) for loc, _, _ in entries)),
            tuple(sorted((word, loc) for loc, word in replaced)),
            tuple(sorted((name, loc) for loc, name, _ in entries)),
        )
        found.append((key, order, replaced, entries))
    best_key = min(key for key, *_ in found)
    _, order, replaced, entries = min(
        (c for c in found if c[0] == best_key), key=lambda c: chain_nuclide_key(mol, c[1], bridges)
    )
    amido = [p for p, w in replaced if w == "amido"]
    if n == 2 and amido and (2, "imido") in replaced and all(p != amido[0] for p, w in replaced if w != "amido"):
        raise UnsupportedStructure("a carbamic acid carrying a carbamoyl group on nitrogen is named as such (P-66.1.6.1.1.4)")
    ends_ligands = [
        w for _, w in replaced if w in _TERMINAL_ONLY
    ]
    if len(ends_ligands) == 2:
        raise UnsupportedStructure("two non-hydroxy end groups make an acyl halide of the polycarbonic acid")
    if numbered is not None:
        numbered.update({center: 2 * i + 1 for i, center in enumerate(order)})
    name = _compose(mol, n, replaced, _letters(mol, graph, order, bridges), _substituent_prefix(entries))
    return insert_chain_nuclides(name, chain_nuclides_text(mol, order, bridges))


def chain_numbering(acid):
    """{carbon: locant} of the numbering that the name of the polycarbonic acid `acid` takes, its nuclides deciding between
    equivalent directions, else {}."""
    from ._isotope_labels import split_isotopes

    split = split_isotopes(acid)
    clean, labels = (split[0], split[1]) if split else (acid, {})
    token = ISOTOPE_LABELS.set({"labels": labels, "consumed": set(), "mol": clean} if labels else None)
    numbered = {}
    try:
        name_polycarbonic(clean, numbered)
    except UnsupportedStructure:
        return {}
    finally:
        ISOTOPE_LABELS.reset(token)
    return numbered


def _locant_number(locant):
    return int(locant.lstrip("N'") or 0)


def _substituent_prefix(entries):
    grouped = {}
    for locant, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    return format_substituent_prefixes(grouped) if grouped else ""


def _acid_group_in(mol, arm):
    """Class 7a acids outrank the polycarbonic acids of class 7b (P-41), so an arm holding one would be the parent."""
    for idx in arm:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetFormalCharge() or atom.GetNumRadicalElectrons():
            return True
        if atom.GetAtomicNum() not in _SYMBOL or not atom.GetTotalNumHs():
            continue
        for x in atom.GetNeighbors():
            if atom.GetAtomicNum() == x.GetAtomicNum() == 8 and x.GetDegree() == 2:
                x = next(m for m in x.GetNeighbors() if m.GetIdx() != atom.GetIdx())
            if x.GetAtomicNum() in (5, 15, 16, 33, 34, 51, 52):
                return True
            if any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(x).GetAtomicNum() in (7, 8, 16, 34, 52) for b in x.GetBonds()):
                return True
    return False


def _arm_names(mol, graph, nitrogen, skip, locant):
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = []
    for root in graph[nitrogen]:
        if root in skip:
            continue
        arm = _arm_atoms(graph, root, nitrogen)
        if arm & skip or _acid_group_in(mol, arm):
            raise UnsupportedStructure("a substituent on a nitrogen of a polycarbonic acid holds a senior group or closes a ring")
        name, compound = name_branch(graph, root, nitrogen, halogens, aromatic, mol=mol, unsaturated=True)
        entries.append((locant, name, compound))
    return entries


def _nitrogen_substituents(mol, graph, order, bridges):
    """P-65.2.3.1.2.1, P-65.2.3.1.4: [(locant, substituent name, is compound)] of the groups on the imido nitrogens
    (N2, N3, ... after the chain position of the nitrogen) and on the amido or hydrazido nitrogen (N, N')."""
    entries = []
    for i, center in enumerate(order):
        word, atom = _oxo_word(mol, center, substituted=True)
        if word == "NH":
            entries += _arm_names(mol, graph, atom, {center}, f"N{2 * i + 1}")
    for i in range(len(order) - 1):
        bridge = bridges[frozenset((order[i], order[i + 1]))]
        if isinstance(bridge, int) and mol.GetAtomWithIdx(bridge).GetAtomicNum() == 7:
            entries += _arm_names(mol, graph, bridge, {order[i], order[i + 1]}, f"N{2 * i + 2}")
    bridge_atoms = _bridge_atoms(bridges)
    for center in (order[0], order[-1]):
        for m in graph[center]:
            atom = mol.GetAtomWithIdx(m)
            if m in bridge_atoms or atom.GetAtomicNum() != 7 or _bond(mol, center, m) != 1.0:
                continue
            betas = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != center]
            entries += _arm_names(mol, graph, m, {center, *betas}, "N")
            for beta in betas:
                entries += _arm_names(mol, graph, beta, {m}, "N'")
    if entries and (
        any(a.GetIsotope() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())
        or any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds())
    ):
        raise UnsupportedStructure("isotopes and stereodescriptors on a substituted polycarbonic acid are not named yet")
    return entries


def _acyl_end_name(mol, graph, order, bridges):
    """P-67.3.1: a polycarbonic acid with an acyloxy end group is formic acid substituted by the rest of the chain,
    '{[(acetyloxy)carbonyl]oxy}formic acid'; None for any other end groups."""
    bridge_atoms = _bridge_atoms(bridges)
    if any(not isinstance(b, int) or mol.GetAtomWithIdx(b).GetAtomicNum() != 8 for b in bridges.values()):
        return None
    if any(_oxo_word(mol, center)[0] != "O" for center in order):
        return None
    ligands = []
    for center in (order[0], order[-1]):
        ligand = [m for m in graph[center] if m not in bridge_atoms and _bond(mol, center, m) == 1.0]
        if len(ligand) != 1:
            return None
        ligands.append(mol.GetAtomWithIdx(ligand[0]))
    hydroxy = [a for a in ligands if a.GetAtomicNum() == 8 and a.GetDegree() == 1 and a.GetTotalNumHs() == 1]
    acyloxy = [a for a in ligands if a.GetAtomicNum() == 8 and a.GetDegree() == 2 and a is not None]
    if len(hydroxy) != 1 or len(acyloxy) != 1:
        return None
    oxygen = acyloxy[0]
    center = order[-1] if ligands[1] is oxygen else order[0]
    (acyl,) = [n.GetIdx() for n in oxygen.GetNeighbors() if n.GetIdx() != center]
    if not any(
        b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(mol.GetAtomWithIdx(acyl)).GetAtomicNum() in (7, 8, 16, 34, 52)
        for b in mol.GetAtomWithIdx(acyl).GetBonds()
    ):
        return None
    chain = order if ligands[0] in hydroxy else order[::-1]
    name, compound = name_branch(graph, acyl, oxygen.GetIdx(), halogen_substituents(mol), mol=mol)
    text, _ = alkoxy_prefix(name, compound)
    for _ in chain[1:]:
        text = enclose(enclose(text) + "carbonyl") + "oxy"
    return enclose(text) + "formic acid"


def _bridge_atoms(bridges):
    atoms = set()
    for bridge in bridges.values():
        atoms.update(bridge if isinstance(bridge, tuple) else (bridge,))
    return atoms


def _replacements(mol, graph, order, bridges, substituted=False):
    """[(locant, prefix word)] and letters for the replaced atoms along the numbered chain."""
    n = len(order)
    replaced, letters = [], []
    for i, center in enumerate(order):
        position = 2 * i + 1
        oxo, _ = _oxo_word(mol, center, substituted)
        if oxo is None:
            raise UnsupportedStructure("unsupported =X atom on a polycarbonic acid")
        if oxo != "O":
            replaced.append((position, _BRIDGE_WORD.get(oxo) or {"NNH2": "hydrazono"}.get(oxo, _INFIX.get(oxo))))
        if i in (0, n - 1):
            bridge_atoms = _bridge_atoms(bridges)
            ligand = [m for m in graph[center] if m not in bridge_atoms and _bond(mol, center, m) == 1.0]
            if len(ligand) != 1:
                raise UnsupportedStructure("a terminal group of the polycarbonic acid is not recognized")
            _terminal(mol, graph, center, ligand[0], position, replaced, letters, substituted)
    for i in range(n - 1):
        bridge = bridges[frozenset((order[i], order[i + 1]))]
        if isinstance(bridge, tuple):
            replaced.append((2 * i + 2, "peroxy"))
            continue
        atom = mol.GetAtomWithIdx(bridge)
        symbol = "NH" if atom.GetAtomicNum() == 7 else _SYMBOL[atom.GetAtomicNum()]
        if symbol != "O":
            replaced.append((2 * i + 2, _BRIDGE_WORD[symbol]))
    return replaced


def _terminal(mol, graph, center, atom_idx, position, replaced, letters, substituted=False):
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
    if atom.GetAtomicNum() in _HALIDE_INFIX and is_halide_ligand(atom):
        replaced.append((position, _PREFIX_WORDS[_HALIDE_INFIX[atom.GetAtomicNum()]]))
        return
    pseudo = pseudohalide_at(mol, atom_idx, center)
    if pseudo is not None and pseudo != "CN":
        replaced.append((position, _PREFIX_WORDS[PSEUDO_INFIX[pseudo]]))
        return
    if atom.GetAtomicNum() == 7 and not atom.GetFormalCharge():
        others = [n for n in atom.GetNeighbors() if n.GetIdx() != center]
        if not others and atom.GetTotalNumHs() == 2:
            replaced.append((position, "amido"))
            return
        if (
            len(others) == 1
            and atom.GetTotalNumHs() == 1
            and others[0].GetAtomicNum() == 7
            and others[0].GetDegree() == 1
            and others[0].GetTotalNumHs() == 2
        ):
            replaced.append((position, "hydrazido"))
            return
        if substituted and not atom.IsInRing() and all(b.GetBondTypeAsDouble() == 1.0 for b in atom.GetBonds()):
            betas = [n for n in others if n.GetAtomicNum() == 7]
            if not betas:
                replaced.append((position, "amido"))
                return
            if len(betas) == 1 and not betas[0].IsInRing() and not betas[0].GetFormalCharge():
                replaced.append((position, "hydrazido"))
                return
    raise UnsupportedStructure("this terminal group of a polycarbonic acid is not supported")


_TERMINAL_ONLY = {"amido", "hydrazido", "fluoro", "chloro", "bromo", "iodo", "azido", "isocyano", "isocyanato", "isothiocyanato", "isoselenocyanato", "isotellurocyanato"}


def _letters(mol, graph, order, bridges):
    """Italic letters with locants of the hydrogen-bearing atoms of the two end groups when a replaced
    chalcogen could sit on either the =X or the -Y position of that carbon (P-65.2.3.1.2.2)."""
    n = len(order)
    found = []
    bridge_atoms = _bridge_atoms(bridges)
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


def _compose(mol, n, replaced, letters="", prefix=""):
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
    separator = "-" if prefix and parts[0][0].isdigit() else ""
    return prefix + separator + "-".join(parts) + core[:-len(" acid")] + (f" {letters}-acid" if letters else " acid")
