"""Substituted nucleosides on the 7 retained names (P-105.2, Chapter P-10).

P-105.2.1: ring substituents (unprimed locants; 'N2'/'O6' on exocyclic atoms),
'thio'/'seleno'/'telluro' replacement, and ribofuranosyl modifications as in
P-102.5 (primed locants): deoxy, deoxy-halo/amino, O-/S-substitution, C-substitution
(P-102.5.6.3, CIP descriptors where required), ring-O thio, O-acyl esters as 'ate'
words (P-102.5.6.1.1), cyclic carbonates (P-101.7.4, P-105.2.2).
P-105.2.2: a principal group senior to the pseudoketone moves the base into a
substituent group or parent ring named by `_ring_group_name`.
Templates are locant-numbered and matched stereo-aware; each sugar X atom may be
absent (deoxy), O/S/N/halogen, or carbon (single or ylidene bond).
"""

import itertools
import re

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._nucleoside import _NUCLEOSIDE_SMILES
from ._numerals import numerical_term
from ._substituents import name_branch, wrap_marks

_RIBO = "[C@H:11]3[C@@H:12]([C@@H:13]([C@H:14]([O:20]3)[CH2:15][OH:25])[OH:23])[OH:22]"
_DEOXY2 = "[C@H:11]2[CH2:12][C@@H:13]([C@H:14]([O:20]2)[CH2:15][OH:25])[OH:23]"
_PARENTS = (
    (
        "adenosine",
        "[CH:2]1=[N:1][C:6](=[C:5]2[C:4](=[N:3]1)[N:9]([CH:8]=[N:7]2)"
        + _RIBO
        + ")[NH2:106]",
        frozenset(),
    ),
    (
        "guanosine",
        "[CH:8]1=[N:7][C:5]2=[C:4]([N:9]1"
        + _RIBO
        + ")[NH:3][C:2](=[N:1][C:6]2=[O:106])[NH2:102]",
        frozenset(),
    ),
    (
        "inosine",
        "[CH:2]1=[N:1][C:6](=[O:106])[C:5]2=[C:4]([NH:3]1)[N:9]([CH:8]=[N:7]2)" + _RIBO,
        frozenset(),
    ),
    (
        "xanthosine",
        "[CH:8]1=[N:7][C:5]2=[C:4]([N:9]1"
        + _RIBO
        + ")[NH:3][C:2](=[O:102])[NH:1][C:6]2=[O:106]",
        frozenset(),
    ),
    (
        "cytidine",
        "[CH:5]1=[CH:6][N:1]([C:2](=[O:102])[N:3]=[C:4]1[NH2:104])" + _RIBO,
        frozenset(),
    ),
    (
        "thymidine",
        "[CH3][C:5]1=[CH:6][N:1]([C:2](=[O:102])[NH:3][C:4]1=[O:104])" + _DEOXY2,
        frozenset({22}),
    ),
    (
        "uridine",
        "[CH:5]1=[CH:6][N:1]([C:2](=[O:102])[NH:3][C:4]1=[O:104])" + _RIBO,
        frozenset(),
    ),
)

_SUGAR_X = {22: ("2", 12), 23: ("3", 13), 25: ("5", 15)}
_STATES = {
    22: ("single", "absent", "double", "gem"),
    23: ("single", "absent", "double", "gem"),
    25: ("single", "absent", "double"),
}
_SUGAR_RING = {11: "1", 12: "2", 13: "3", 14: "4", 15: "5"}
_EXOCYCLIC = {102: "2", 104: "4", 106: "6"}
_BASE_RING = frozenset(range(1, 10))
_CHALCOGEN = {16: "thio", 34: "seleno", 52: "telluro"}
_OXO_PREFIX = {8: "oxo", 16: "sulfanylidene", 34: "selanylidene", 52: "tellanylidene"}
_ANION_ALIASES = {"ethanoate": "acetate"}
_ALLOWED_Z = {6, 7, 8, 9, 16, 17, 34, 35, 52, 53}
_SCAFFOLD = Chem.MolFromSmarts("[n]C1CCC[O,S,Se,Te]1")
_BIS = {2: "bis", 3: "tris", 4: "tetrakis", 5: "pentakis", 6: "hexakis"}


class _Entry:
    __slots__ = ("compound", "hetero", "key", "locant", "name", "position", "prime")

    def __init__(self, locant, prime, name, compound=False, hetero="", key=None):
        self.locant = locant
        self.prime = prime
        self.name = name
        self.compound = compound
        self.hetero = hetero
        self.key = key
        self.position = int(re.search(r"\d+", locant).group())

    def text(self, primed=None):
        return self.locant + ("′" if (self.prime if primed is None else primed) else "")


def _strip_maps(mol):
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return copy


def _replace(rw, idx, smarts, map_num):
    query_atom = Chem.AtomFromSmarts(smarts)
    query_atom.SetAtomMapNum(map_num)
    rw.ReplaceAtom(idx, query_atom, preserveProps=True)


def _make_query(parent, plan):
    rw = Chem.RWMol(parent)
    by_map = {a.GetAtomMapNum(): a.GetIdx() for a in rw.GetAtoms() if a.GetAtomMapNum()}
    _replace(rw, by_map[20], "[O,S,Se,Te]", 20)
    for map_num in _EXOCYCLIC:
        if map_num in by_map and rw.GetAtomWithIdx(by_map[map_num]).GetAtomicNum() == 8:
            _replace(rw, by_map[map_num], "[O,S,Se,Te]", map_num)
    removed = []
    for map_num, state in plan.items():
        x_idx = by_map[map_num]
        carbon = rw.GetAtomWithIdx(x_idx).GetNeighbors()[0]
        if state in ("single", "gem"):
            _replace(rw, x_idx, "[O,S,N,F,Cl,Br,I,#6]", map_num)
            if state == "gem":
                carbon.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            continue
        carbon.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if state == "double":
            _replace(rw, x_idx, "[#6]", map_num)
            rw.GetBondBetweenAtoms(carbon.GetIdx(), x_idx).SetBondType(
                Chem.BondType.DOUBLE
            )
        else:
            removed.append(x_idx)
    for x_idx in sorted(removed, reverse=True):
        rw.RemoveAtom(x_idx)
    return rw.GetMol()


def _build_queries():
    queries = []
    for order, (name, mapped, fixed_absent) in enumerate(_PARENTS):
        parent = Chem.MolFromSmiles(mapped)
        assert Chem.MolToSmiles(_strip_maps(parent)) == Chem.CanonSmiles(
            _NUCLEOSIDE_SMILES[name]
        ), name
        variable = [m for m in _SUGAR_X if m not in fixed_absent]
        for states in itertools.product(*[_STATES[m] for m in variable]):
            plan = dict(zip(variable, states))
            queries.append((order, name, fixed_absent, plan, _make_query(parent, plan)))
    return queries


_QUERIES = _build_queries()


def _branch_atoms(graph, root, came_from):
    seen = {root}
    stack = [root]
    while stack:
        for nb in graph[stack.pop()]:
            if nb != came_from and nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return seen


def _name_substituent(
    mol, graph, halogens, aromatic, core, root, came_from, double_ok=False
):
    if (
        mol.GetBondBetweenAtoms(root, came_from).GetBondTypeAsDouble() != 1.0
        and not double_ok
    ):
        raise UnsupportedStructure(
            "non-single bond from a nucleoside core atom to its substituent"
        )
    atoms = _branch_atoms(graph, root, came_from)
    if atoms & core:
        raise UnsupportedStructure("substituent reconnects to the nucleoside core")
    if any(
        mol.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        for i in atoms
    ):
        raise UnsupportedStructure("stereo substituent on a nucleoside")
    for b in mol.GetBonds():
        if (
            b.GetBeginAtomIdx() in atoms
            and b.GetEndAtomIdx() in atoms
            and b.GetStereo() != Chem.BondStereo.STEREONONE
        ):
            raise UnsupportedStructure("stereo substituent on a nucleoside")
    if _is_azide(mol, atoms, root):
        return "azido", False
    name, compound = name_branch(graph, root, came_from, halogens, aromatic, mol)
    return name.replace("ethanoyl", "acetyl"), compound


def _is_azide(mol, atoms, root):
    if len(atoms) != 3 or sorted(
        mol.GetAtomWithIdx(i).GetAtomicNum() for i in atoms
    ) != [7, 7, 7]:
        return False
    return sorted(mol.GetAtomWithIdx(i).GetFormalCharge() for i in atoms) == [-1, 0, 1]


def _cip_data(mol):
    ranked = Chem.Mol(mol)
    Chem.AssignStereochemistry(ranked, cleanIt=True, force=True)
    ranks = {
        a.GetIdx(): a.GetIntProp("_CIPRank") if a.HasProp("_CIPRank") else 0
        for a in ranked.GetAtoms()
    }
    labelled = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(labelled)
    atoms = {
        a.GetIdx(): a.GetProp("_CIPCode")
        for a in labelled.GetAtoms()
        if a.HasProp("_CIPCode")
    }
    bonds = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())): b.GetProp("_CIPCode")
        for b in labelled.GetBonds()
        if b.HasProp("_CIPCode")
    }
    return ranks, atoms, bonds


def _acid_anion(mol, graph, o_idx, acyl_root):
    from .core import _name_mol

    atoms = _branch_atoms(graph, acyl_root, o_idx) | {o_idx}
    if any(
        mol.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        for i in atoms
    ):
        raise UnsupportedStructure("stereo acyl group on a nucleoside")
    acid = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, sorted(atoms)))
    name = _name_mol(acid) if acid is not None else ""
    if not name.endswith("ic acid"):
        raise UnsupportedStructure("acyl group is not a nameable carboxylic acid")
    return _ANION_ALIASES.get(
        name[: -len("ic acid")] + "ate", name[: -len("ic acid")] + "ate"
    )


def _is_carboxylic_acyl(mol, root):
    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
        return False
    bonds = [
        (
            n.GetAtomicNum(),
            mol.GetBondBetweenAtoms(root, n.GetIdx()).GetBondTypeAsDouble(),
        )
        for n in atom.GetNeighbors()
    ]
    return (8, 2.0) in bonds and sum(
        1 for z, order in bonds if z == 6 and order == 1.0
    ) == 1


def _is_plain_carbonate_carbon(mol, root):
    atom = mol.GetAtomWithIdx(root)
    return atom.GetAtomicNum() == 6 and sorted(
        n.GetAtomicNum() for n in atom.GetNeighbors()
    ) == [8, 8, 8]


class _Analysis:
    def __init__(self, parent, order, absent, by_map, core):
        self.parent = parent
        self.order = order
        self.absent = absent
        self.by_map = by_map
        self.core = core
        self.entries = []
        self.esters = []
        self.descriptors = []
        self.carbonate = None

    @property
    def cost(self):
        return len(self.entries) + len(self.esters) + len(self.descriptors)


def _sanitized(mol, kekulize=False):
    try:
        if kekulize:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        else:
            Chem.SanitizeMol(mol)
    except ValueError:
        return None
    return mol


def _lactim_ethers(mol):
    """Yield (mol2, {O idx: R idx}) with each O-substituted lactim oxygen
    rewritten as a lactam C=O so the oxo templates match."""
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 2 or atom.IsInRing():
            continue
        carbons = [
            n
            for n in atom.GetNeighbors()
            if n.GetIsAromatic() and n.GetAtomicNum() == 6
        ]
        others = [n for n in atom.GetNeighbors() if not n.GetIsAromatic()]
        if len(carbons) != 1 or len(others) != 1 or others[0].GetAtomicNum() != 6:
            continue
        nitrogens = [
            n.GetIdx()
            for n in carbons[0].GetNeighbors()
            if n.GetIsAromatic()
            and n.GetAtomicNum() == 7
            and n.GetDegree() == 2
            and n.GetTotalNumHs() == 0
        ]
        if nitrogens:
            candidates.append(
                (atom.GetIdx(), carbons[0].GetIdx(), others[0].GetIdx(), nitrogens)
            )
    if not candidates:
        return
    for choice in itertools.product(*[c[3] for c in candidates]):
        if len(set(choice)) != len(choice):
            continue
        rw = Chem.RWMol(mol)
        ethers = {}
        for (o, c, r, _), n in zip(candidates, choice):
            rw.RemoveBond(o, r)
            rw.GetBondBetweenAtoms(o, c).SetBondType(Chem.BondType.DOUBLE)
            rw.GetAtomWithIdx(n).SetNumExplicitHs(1)
            rw.GetAtomWithIdx(n).SetNoImplicit(True)
            ethers[o] = r
        out = _sanitized(rw.GetMol())
        if out is not None:
            yield out, ethers


def _imine_variants(mol):
    """Yield (mol2, {ring N idx: R idx}) with each N-substituted ring nitrogen next
    to an exocyclic C=NH rewritten as the unsubstituted amino-azine tautomer."""
    candidates = []
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() != 7
            or not atom.IsInRing()
            or atom.GetDegree() != 3
            or atom.GetTotalNumHs()
            or atom.GetFormalCharge()
        ):
            continue
        subs = [
            n
            for n in atom.GetNeighbors()
            if not mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).IsInRing()
        ]
        if len(subs) != 1 or subs[0].GetAtomicNum() != 6:
            continue
        for c in atom.GetNeighbors():
            if (
                not mol.GetBondBetweenAtoms(atom.GetIdx(), c.GetIdx()).IsInRing()
                or c.GetAtomicNum() != 6
            ):
                continue
            for e in c.GetNeighbors():
                if (
                    e.GetAtomicNum() == 7
                    and e.GetDegree() == 1
                    and mol.GetBondBetweenAtoms(
                        c.GetIdx(), e.GetIdx()
                    ).GetBondTypeAsDouble()
                    == 2.0
                ):
                    candidates.append(
                        (atom.GetIdx(), subs[0].GetIdx(), c.GetIdx(), e.GetIdx())
                    )
    for n, r, c, e in candidates:
        rw = Chem.RWMol(mol)
        if _sanitized(rw, kekulize=True) is None:
            continue
        rw.RemoveBond(n, r)
        rw.GetBondBetweenAtoms(c, e).SetBondType(Chem.BondType.SINGLE)
        rw.GetBondBetweenAtoms(n, c).SetBondType(Chem.BondType.DOUBLE)
        out = _sanitized(rw.GetMol())
        if out is not None:
            yield out, {n: r}


def _evaluate(ctx, match, query, parent, order, fixed_absent, plan, ethers, ring_subs):
    mol, graph, graph2, halogens, aromatic, (ranks, atom_cip, bond_cip) = ctx
    by_map = {
        a.GetAtomMapNum(): match[a.GetIdx()]
        for a in query.GetAtoms()
        if a.GetAtomMapNum()
    }
    core = set(match)
    absent = frozenset(fixed_absent) | {m for m, s in plan.items() if s == "absent"}
    an = _Analysis(parent, order, absent, by_map, core)
    for idx in core:
        atom = mol.GetAtomWithIdx(idx)
        if (
            atom.GetAtomicNum() not in _ALLOWED_Z
            or atom.GetFormalCharge()
            or atom.GetIsotope()
            or atom.GetNumRadicalElectrons()
        ):
            raise UnsupportedStructure(
                "charged, isotopic or unsupported atom in a nucleoside"
            )
    extras = lambda idx: [n for n in graph2[idx] if n not in core]
    for query_atom in query.GetAtoms():
        if not query_atom.GetAtomMapNum() and extras(match[query_atom.GetIdx()]):
            raise UnsupportedStructure("substituent on an unlabelled core atom")

    def add(root, came_from, locant, prime, hetero="", scope=None, double_ok=False):
        name, compound = _name_substituent(
            mol, graph, halogens, aromatic, scope or core, root, came_from, double_ok
        )
        an.entries.append(_Entry(locant, prime, name, compound, hetero))

    def descriptor(carbon_idx, c_map):
        label = atom_cip.get(carbon_idx)
        if label:
            an.descriptors.append((int(_SUGAR_RING[c_map]), label))

    for map_num in sorted(by_map):
        idx = by_map[map_num]
        z = mol.GetAtomWithIdx(idx).GetAtomicNum()
        if map_num in _BASE_RING:
            for nb in extras(idx):
                add(nb, idx, str(map_num), False)
            if idx in ring_subs:
                add(ring_subs[idx], idx, str(map_num), False)
        elif map_num in _EXOCYCLIC:
            position = _EXOCYCLIC[map_num]
            if z == 7:
                for nb in extras(idx):
                    add(nb, idx, f"N{position}", False)
            elif extras(idx):
                raise UnsupportedStructure("substituted oxo group on a nucleoside base")
            elif z in _CHALCOGEN:
                an.entries.append(_Entry(position, False, _CHALCOGEN[z]))
            elif idx in ethers:
                add(ethers[idx], idx, f"O{position}", False)
        elif map_num == 20 and z in _CHALCOGEN:
            an.entries.append(_Entry("4", True, _CHALCOGEN[z]))
        elif map_num in (11, 14) and extras(idx):
            for nb in extras(idx):
                add(nb, idx, _SUGAR_RING[map_num], True, "C")
            descriptor(idx, map_num)

    carbonyls = {}
    for x_map, (position, c_map) in _SUGAR_X.items():
        c, x, state = by_map[c_map], by_map.get(x_map), plan.get(x_map)
        sigma = extras(c)
        if x is None:
            if sigma:
                raise UnsupportedStructure("replacement group on a deoxy position")
            if x_map in absent and x_map not in fixed_absent:
                an.entries.append(_Entry(position, True, "deoxy"))
            continue
        if state == "gem":
            if len(sigma) != 1 or ranks[x] != ranks[sigma[0]]:
                raise UnsupportedStructure(
                    "a free sugar carbon must carry two identical groups"
                )
        elif any(ranks[x] <= ranks[e] for e in sigma):
            raise UnsupportedStructure(
                "the senior group at a sugar carbon must stand for the hydroxy group"
            )
        z = mol.GetAtomWithIdx(x).GetAtomicNum()
        if mol.GetBondBetweenAtoms(c, x).GetBondTypeAsDouble() == 2.0:
            an.entries.append(_Entry(position, True, "deoxy"))
            add(x, c, position, True, scope=core - {x}, double_ok=True)
            label = bond_cip.get(frozenset((c, x)))
            if label:
                an.descriptors.append((int(position), label))
            continue
        if z in (8, 16):
            if z == 16:
                an.entries.append(_Entry(position, True, "thio"))
            for nb in extras(x):
                if z == 8 and (
                    _is_carboxylic_acyl(mol, nb) or _is_plain_carbonate_carbon(mol, nb)
                ):
                    carbonyls.setdefault(nb, []).append((position, x))
                else:
                    add(nb, x, position, True, "O" if z == 8 else "S", core | {x})
            for nb in sigma:
                add(nb, c, position, True, "C")
        else:
            an.entries.append(_Entry(position, True, "deoxy"))
            if z in (9, 17, 35, 53):
                an.entries.append(_Entry(position, True, halogens[x]))
            elif z == 7 and not extras(x):
                an.entries.append(_Entry(position, True, "amino"))
            else:
                add(x, c, position, True, scope=core - {x})
            for nb in sigma:
                add(nb, c, position, True)
        if sigma and state != "gem" and (z not in (8, 16) or x_map == 25):
            descriptor(c, c_map)

    for root, users in carbonyls.items():
        if _is_plain_carbonate_carbon(mol, root):
            if len(users) != 2 or len(carbonyls) != 1:
                raise UnsupportedStructure(
                    "unsupported carbonate ester on a nucleoside"
                )
            an.carbonate = sorted(int(p) for p, _ in users)
            an.entries.extend(_Entry(p, True, "deoxy") for p, _ in users)
        elif len(users) == 1:
            position, x = users[0]
            an.esters.append((int(position), _acid_anion(mol, graph, x, root), root, x))
        else:
            raise UnsupportedStructure(
                "acyl bridge between two nucleoside hydroxy groups"
            )
    if an.carbonate and an.esters:
        raise UnsupportedStructure("carbonate together with esters on a nucleoside")
    return an


def _alpha_key(entry):
    if entry.key:
        return entry.key
    return re.sub(r"[^a-z]", "", re.sub(r"^[\d,′\-]+", "", entry.name.lower()))


def _group_text(entries, primed=None):
    groups = {}
    for entry in entries:
        groups.setdefault(
            (entry.hetero, entry.name, entry.compound, _alpha_key(entry)), []
        ).append(entry)
    parts = []
    for (hetero, name, compound, key), members in groups.items():
        members.sort(key=lambda e: (e.position, e.locant))
        count = len(members)
        multiplier = (
            (_BIS[count] if compound else numerical_term(count)) if count > 1 else ""
        )
        label = f"{hetero}-" if hetero else ""
        joiner = "-" if multiplier and hetero else ""
        locants = ",".join(m.text(primed) for m in members)
        parts.append(
            (
                key,
                f"{locants}-{multiplier}{joiner}{label}{wrap_marks(name) if compound else name}",
            )
        )
    return [text for _, text in sorted(parts)]


def _ester_words(esters):
    by_anion = {}
    for position, anion, *_ in esters:
        by_anion.setdefault(anion, []).append(position)
    words = []
    for anion in sorted(by_anion):
        count = len(by_anion[anion])
        compound = bool(re.search(r"[\d ]", anion))
        multiplier = (
            (_BIS[count] if compound else numerical_term(count)) if count > 1 else ""
        )
        shown = wrap_marks(anion) if compound and count > 1 else anion
        words.append(
            f"{','.join(f'{p}′' for p in sorted(by_anion[anion]))}-{multiplier}{shown}"
        )
    return " ".join(words)


def _plain_name(an):
    labels = ",".join(f"{pos}′{label}" for pos, label in sorted(set(an.descriptors)))
    name = (
        (f"({labels})-" if labels else "")
        + "-".join(_group_text(an.entries))
        + an.parent
    )
    if an.carbonate:
        name += f"-{an.carbonate[0]}′,{an.carbonate[1]}′-diyl carbonate"
    if an.esters:
        name += " " + _ester_words(an.esters)
    return name


_PYRIMIDINE_BONDS = ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1))
_PURINE_BONDS = _PYRIMIDINE_BONDS + ((4, 9), (5, 7), (7, 8), (8, 9))
_RING_STEM = {
    "purine": ("purin", _PURINE_BONDS),
    "pyrimidine": ("pyrimidin", _PYRIMIDINE_BONDS),
}
_SWAP = {1: 3, 3: 1, 4: 6, 6: 4}
_SUFFIXES = {
    "acid": "carboxylic acid",
    "ester": "carboxylate",
    "amide": "carboxamide",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
}
_SENIORITY = tuple(_SUFFIXES)


def _has_perfect_matching(nodes, bonds):
    if not nodes:
        return True
    first = min(nodes)
    return any(
        first in (a, b) and _has_perfect_matching(nodes - {a, b}, bonds)
        for a, b in bonds
        if a in nodes and b in nodes
    )


def _ring_group_name(kind, prefix_entries, saturated, free_valence=None, suffix=None):
    """Name the purine/pyrimidine base as a substituent group ('...-5-yl') or,
    with `suffix=(positions, class)`, as the parent ring of a senior group.
    Numbering follows P-31.1.4: indicated hydrogen, free valence/suffix, hydro,
    all detachable prefixes, first-cited prefix; hydro prefixes cover the
    saturated positions (P-58.2.3, P-58.2.4)."""
    stem, bonds = _RING_STEM[kind]
    nodes = {n for b in bonds for n in b}
    valid = (
        {p for p in nodes if _has_perfect_matching(nodes - {p}, bonds)}
        if kind == "purine"
        else set()
    )
    best = None
    for swapped in (False, True) if kind == "pyrimidine" else (False,):
        remap = (lambda p: _SWAP.get(p, p)) if swapped else (lambda p: p)
        sat = sorted(remap(p) for p in saturated)
        indicated = min((p for p in sat if p in valid), default=None)
        if kind == "purine" and indicated is None:
            continue
        hydro = [p for p in sat if p != indicated]
        if len(hydro) % 2:
            continue
        entries = [(remap(p), e) for p, e in prefix_entries]
        marker = sorted(remap(p) for p in (suffix[0] if suffix else [free_valence]))
        first = min(entries, key=lambda pe: _alpha_key(pe[1]))[0] if entries else 0
        key = (
            indicated or 0,
            tuple(marker),
            tuple(hydro),
            tuple(sorted(p for p, _ in entries)),
            first,
        )
        if best is None or key < best[0]:
            best = (key, entries, indicated, hydro, marker)
    if best is None:
        raise UnsupportedStructure(
            "no consistent hydro/indicated-hydrogen numbering for the base"
        )
    _, entries, indicated, hydro, marker = best
    grouped = {}
    for p, e in entries:
        grouped.setdefault((e.name, e.compound, _alpha_key(e)), []).append(p)
    parts = []
    for (name, compound, key), positions in grouped.items():
        count = len(positions)
        multiplier = (
            (_BIS[count] if compound else numerical_term(count)) if count > 1 else ""
        )
        parts.append(
            (
                key,
                f"{','.join(map(str, sorted(positions)))}-{multiplier}{wrap_marks(name) if compound else name}",
            )
        )
    core = (
        f"{','.join(map(str, hydro))}-{numerical_term(len(hydro))}hydro"
        if hydro
        else ""
    )
    if indicated:
        core += ("-" if hydro else "") + f"{indicated}H-"
    text = "".join(f"{t}-" for _, t in sorted(parts)) + core
    locants = ",".join(map(str, marker))
    if suffix:
        multiplier = numerical_term(len(marker)) if len(marker) > 1 else ""
        return f"{text}{stem}e-{locants}-{multiplier}{_SUFFIXES[suffix[1]]}"
    return f"{text}{stem}-{locants}-yl"


def _glycosyl_group(an, name_acyl):
    deoxy_2, deoxy_3 = 22 in an.absent, 23 in an.absent
    stem = (
        "β-D-glycero-pentofuranosyl"
        if deoxy_2 and deoxy_3
        else "β-D-erythro-pentofuranosyl"
        if deoxy_2 or deoxy_3
        else "β-D-ribofuranosyl"
    )
    if an.carbonate or an.descriptors:
        raise UnsupportedStructure(
            "unsupported sugar modification on a substituent nucleoside group"
        )
    entries = [e for e in an.entries if e.prime]
    for position, _, root, x in an.esters:
        name, compound = name_acyl(root, x)
        entries.append(_Entry(str(position), True, name, compound, "O"))
    prefixes = _group_text(entries, primed=False)
    text = "-".join(prefixes) + ("-" if prefixes else "") + stem
    key = re.sub(
        r"[^a-z]", "", re.sub(r"^[\d,\-]+", "", re.sub(r"[αβ]-d-", "", text.lower()))
    )
    return text, bool(prefixes), key


_SENIOR_SMARTS = (
    ("acid", Chem.MolFromSmarts("[CX3](=O)[OX2H1]")),
    ("ester", Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")),
    ("amide", Chem.MolFromSmarts("[CX3](=O)[NX3]")),
    ("nitrile", Chem.MolFromSmarts("[CX2]#[NX1]")),
    ("aldehyde", Chem.MolFromSmarts("[CX3H1]=O")),
    ("sulfonic acid", Chem.MolFromSmarts("[SX4](=O)(=O)[OX2H1]")),
    ("acyl halide", Chem.MolFromSmarts("[CX3](=O)[F,Cl,Br,I]")),
)


def _seniors_in(mol, atoms, root, came_from):
    """Senior groups whose carbon lies in `atoms`; an acyl group directly on
    an exocyclic N/O (N-acyl) stays an N-/O-substituent."""
    found = []
    for kind, pattern in _SENIOR_SMARTS:
        for match in mol.GetSubstructMatches(pattern):
            carbon = mol.GetAtomWithIdx(match[0])
            if match[0] not in atoms or carbon.IsInRing():
                continue
            root_atom = mol.GetAtomWithIdx(root)
            if (
                match[0] == root
                and mol.GetAtomWithIdx(came_from).GetAtomicNum() in (7, 8)
            ) or (
                root_atom.GetAtomicNum() in (7, 8)
                and any(n.GetIdx() == match[0] for n in root_atom.GetNeighbors())
            ):
                continue
            found.append((kind, match[0]))
    return found


def _senior_name(mol, graph, halogens, aromatic, an):
    from .core import _name_mol

    by_map = an.by_map
    kind = (
        "purine"
        if an.parent in ("adenosine", "guanosine", "inosine", "xanthosine")
        else "pyrimidine"
    )
    ring = {p: by_map[p] for p in _BASE_RING if p in by_map}
    ring_atoms = set(ring.values())
    sugar_atoms = {
        by_map[m] for m in by_map if m in _SUGAR_RING or m in (20, 22, 23, 25)
    }
    core = ring_atoms | sugar_atoms
    for idx in sugar_atoms:
        for nb in graph[idx]:
            if nb not in core and _seniors_in(
                mol, _branch_atoms(graph, nb, idx), nb, idx
            ):
                raise UnsupportedStructure("senior group on a sugar substituent")
    branches = []
    for pos, idx in ring.items():
        for nb in graph[idx]:
            if nb not in ring_atoms and nb not in sugar_atoms:
                atoms = _branch_atoms(graph, nb, idx)
                branches.append((pos, idx, nb, atoms, _seniors_in(mol, atoms, nb, idx)))
    all_found = {k for *_, found in branches for k, _ in found}
    if not all_found:
        return None
    top = next((c for c in _SENIORITY if c in all_found), None)
    if top is None or all_found - set(_SENIORITY):
        raise UnsupportedStructure(
            "senior group class not supported on a nucleoside base"
        )
    sugar_n = next(p for p, idx in ring.items() if by_map[11] in graph[idx])
    saturated = {
        p
        for p, idx in ring.items()
        if (
            mol.GetAtomWithIdx(idx).GetAtomicNum() == 7
            and mol.GetAtomWithIdx(idx).GetDegree()
            + mol.GetAtomWithIdx(idx).GetTotalNumHs()
            == 3
        )
        or any(
            mol.GetBondBetweenAtoms(idx, n).GetBondTypeAsDouble() == 2.0
            and n not in ring_atoms
            for n in graph[idx]
        )
    }
    prefixes, tops = [], []
    for pos, idx, nb, atoms, found in branches:
        atom = mol.GetAtomWithIdx(nb)
        if (
            mol.GetBondBetweenAtoms(idx, nb).GetBondTypeAsDouble() == 2.0
            and atom.GetDegree() == 1
            and atom.GetAtomicNum() in _OXO_PREFIX
        ):
            prefixes.append(
                (pos, _Entry(str(pos), False, _OXO_PREFIX[atom.GetAtomicNum()]))
            )
        elif any(k == top for k, _ in found):
            tops.append((pos, nb, atoms, found))
        else:
            name, compound = _name_substituent(
                mol, graph, halogens, aromatic, core, nb, idx
            )
            prefixes.append((pos, _Entry(str(pos), False, name, compound)))
    glycosyl, compound, key = _glycosyl_group(
        an,
        lambda root, x: _name_substituent(
            mol, graph, halogens, aromatic, core, root, x
        ),
    )
    prefixes.append((sugar_n, _Entry(str(sugar_n), False, glycosyl, compound, key=key)))
    if all(any(k == top and c == nb for k, c in found) for _, nb, _, found in tops):
        alcohol = None
        for pos, nb, _, _ in tops:
            if top == "ester":
                o = next(
                    n
                    for n in graph[nb]
                    if mol.GetAtomWithIdx(n).GetAtomicNum() == 8
                    and mol.GetBondBetweenAtoms(nb, n).GetBondTypeAsDouble() == 1.0
                )
                r = next(n for n in graph[o] if n != nb)
                found_alcohol, _ = _name_substituent(
                    mol, graph, halogens, aromatic, core | {nb, o}, r, o
                )
                if alcohol not in (None, found_alcohol):
                    raise UnsupportedStructure("different ester groups on one base")
                alcohol = found_alcohol
            if top == "amide" and any(
                mol.GetAtomWithIdx(n).GetAtomicNum() == 7
                and mol.GetAtomWithIdx(n).GetDegree() != 1
                for n in graph[nb]
            ):
                raise UnsupportedStructure(
                    "substituted carboxamide on a nucleoside base"
                )
        name = _ring_group_name(
            kind, prefixes, saturated, suffix=([p for p, *_ in tops], top)
        )
        return f"{alcohol} {name}" if alcohol else name
    if len(tops) != 1:
        raise UnsupportedStructure(
            "senior groups in more than one place on a nucleoside base"
        )
    pos, nb, atoms, _ = tops[0]
    group = _ring_group_name(kind, prefixes, saturated, free_valence=pos)
    rw = Chem.RWMol(mol)
    placeholder = rw.AddAtom(Chem.Atom(53))
    rw.AddBond(nb, placeholder, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(placeholder).SetProp("_named_prefix", wrap_marks(group))
    for idx in sorted(set(range(mol.GetNumAtoms())) - atoms, reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return _name_mol(out)


_CACHE = {}


def _best_name(mol):
    key = Chem.MolToSmiles(mol)
    if key not in _CACHE:
        _CACHE[key] = (
            _compute(mol)
            if len(Chem.GetMolFrags(mol)) == 1 and mol.HasSubstructMatch(_SCAFFOLD)
            else None
        )
    return _CACHE[key]


def _compute(mol):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    cip = _cip_data(mol)
    best = None
    variants = [(mol, {}, {})]
    variants += [(m2, e, {}) for m2, e in _lactim_ethers(mol)] + [
        (m2, {}, r) for m2, r in _imine_variants(mol)
    ]
    for mol2, ethers, ring_subs in variants:
        ctx = (mol, graph, adjacency(mol2), halogens, aromatic, cip)
        for order, parent, fixed_absent, plan, query in _QUERIES:
            for match in mol2.GetSubstructMatches(
                query, useChirality=True, uniquify=False, maxMatches=64
            ):
                try:
                    an = _evaluate(
                        ctx,
                        match,
                        query,
                        parent,
                        order,
                        fixed_absent,
                        plan,
                        ethers,
                        ring_subs,
                    )
                    if best is not None and (an.cost, order) > best[0]:
                        continue
                    senior = _senior_name(mol, graph, halogens, aromatic, an)
                    name = _plain_name(an) if senior is None else senior
                except UnsupportedStructure:
                    continue
                if best is None or (an.cost, order, name) < (*best[0], best[1]):
                    best = ((an.cost, order), name)
    return None if best is None else best[1]


def has_substituted_nucleoside_name(mol) -> bool:
    return _best_name(mol) is not None


def name_substituted_nucleoside(mol) -> str:
    name = _best_name(mol)
    if name is None:
        raise UnsupportedStructure("not a substituted nucleoside")
    return name
