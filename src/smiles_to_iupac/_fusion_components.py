"""Fusion components of P-25.1 and P-25.2: the retained and systematic ring systems that can be a parent or an
attached component, recognised from the ring subset of a fused system by their skeleton graph."""

from dataclasses import dataclass, field
from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fused_numbering import HETERO_ORDER, HETERO_RANK, _locant_key, fused_numberings, orientation_key
from ._numerals import numerical_term

PARENT_HETERO_ORDER = ("N", "F", "Cl", "Br", "I", "O", "S", "Se", "Te", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B", "Al", "Ga", "In", "Tl")
PARENT_HETERO_RANK = {e: i for i, e in enumerate(PARENT_HETERO_ORDER)}
A_PREFIX = {
    "O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha", "As": "arsa", "Sb": "stiba",
    "Bi": "bisma", "Si": "sila", "Ge": "germa", "Sn": "stanna", "Pb": "plumba", "B": "bora", "Al": "aluma", "Ga": "galla",
    "In": "indiga", "Tl": "thalla", "F": "fluora", "Cl": "chlora", "Br": "broma", "I": "ioda",
}
_SIX_A = {"O", "S", "Se", "Te", "Bi"}
_SIX_B = {"N", "Si", "Ge", "Sn", "Pb"}
_STEMS = {3: "irene", 4: "ete", 5: "ole", 7: "epine", 8: "ocine", 9: "onine", 10: "ecine"}
_RETAINED_MONOCYCLES = {
    (5, ("O",)): "furan", (5, ("S",)): "thiophene", (5, ("Se",)): "selenophene", (5, ("Te",)): "tellurophene", (5, ("N",)): "pyrrole",
    (5, ("N", "N")): None, (6, ("N",)): "pyridine", (6, ("N", "N")): None, (6, ("O",)): "pyran", (6, ("S",)): "thiopyran",
    (6, ("Se",)): "selenopyran", (6, ("Te",)): "telluropyran",
}
_RETAINED_PREFIX = {
    "furan": "furo", "thiophene": "thieno", "selenophene": "selenopheno", "tellurophene": "telluropheno", "pyrrole": "pyrrolo",
    "imidazole": "imidazo", "pyrazole": "pyrazolo", "pyridine": "pyrido", "pyrimidine": "pyrimido", "pyridazine": "pyridazino",
    "pyrazine": "pyrazino", "pyran": "pyrano", "thiopyran": "thiopyrano", "selenopyran": "selenopyrano", "telluropyran": "telluropyrano",
}


@dataclass
class Component:
    """A fusion component bound to atoms of a ring system: `numberings` are dicts atom -> locant."""

    name: str
    prefix: str
    kind: str
    atoms: list
    numberings: list
    elements: list
    ring_sizes: list
    retained: bool = True
    benzo_unit: bool = False
    hetero_text: str = ""
    senior_key: tuple = ()
    orientation: tuple = ()
    extra: dict = field(default_factory=dict)

    def bracketed_prefix(self):
        return (f"[{self.hetero_text}]" if self.hetero_text else "") + self.prefix

    def full_name(self):
        return (f"{self.hetero_text}-" if self.hetero_text else "") + self.name


@dataclass
class Prototype:
    name: str
    prefix: str
    kind: str
    graph: Chem.Mol
    retained: bool = True
    benzo_unit: bool = False
    hetero_text: str = ""
    fixed_numbering: dict = None
    exception: str = None


def skeleton(atom_symbols, bonds):
    rw = Chem.RWMol()
    for symbol in atom_symbols:
        atom = Chem.Atom(symbol)
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
    for a, b in bonds:
        rw.AddBond(a, b, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    return mol


def skeleton_of(mol):
    return skeleton([a.GetSymbol() for a in mol.GetAtoms()], [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])


def graph_key(mol):
    return Chem.MolToSmiles(mol, canonical=True)


def _from_smiles(smiles):
    return skeleton_of(Chem.MolFromSmiles(smiles, sanitize=False))


def chain_graph(sizes, offsets):
    """Ortho-fused chain of rings; ring k+1 is fused to the edge of ring k that lies offsets[k-1] edges clockwise from
    the edge shared with ring k-1 (3 is straight in a hexagon, 2 is angular)."""
    atoms = ["C"] * sizes[0]
    bonds = [(i, (i + 1) % sizes[0]) for i in range(sizes[0])]
    previous = list(range(sizes[0]))
    edge = 0
    for k in range(1, len(sizes)):
        if k > 1:
            edge = (0 + offsets[k - 2]) % len(previous)
        a, b = previous[edge], previous[(edge + 1) % len(previous)]
        new = list(range(len(atoms), len(atoms) + sizes[k] - 2))
        atoms += ["C"] * len(new)
        cycle = [b, a] + new
        bonds += [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(1, len(cycle))]
        previous = cycle
    return skeleton(atoms, sorted({tuple(sorted(x)) for x in bonds}))


def _central_ring_graph(central, chains):
    atoms = ["C"] * central
    bonds = [(i, (i + 1) % central) for i in range(central)]
    for index, sizes in enumerate(chains):
        shared = (2 * index + 1, 2 * index)
        for size in sizes:
            new = list(range(len(atoms), len(atoms) + size - 2))
            atoms += ["C"] * len(new)
            cycle = [shared[0], shared[1]] + new
            bonds += [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(1, len(cycle))]
            opposite = size // 2
            shared = (cycle[opposite + 1], cycle[opposite])
    return skeleton(atoms, sorted({tuple(sorted(x)) for x in bonds}))


def _elide(stem, ending):
    return stem[:-1] + ending if stem[-1] == ending[0] else stem + ending


_HYDROCARBONS = [
    ("ovalene", "C1=CC2=C3C4=C1C=CC5=CC6=C7C8=C(C=CC9=C8C1=C(C=C9)C=C(C3=C1C7=C54)C=C2)C=C6"),
    ("pyranthrene", "C1=CC=C2C(=C1)C=C3C=CC4=CC5=C6C(=CC7=CC=CC=C75)C=CC8=CC2=C3C4=C86"),
    ("coronene", "C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7"),
    ("rubicene", "C1=CC=C2C(=C1)C3=C4C2=C5C=CC=C6C5=C(C4=CC=C3)C7=CC=CC=C76"),
    ("perylene", "C1=CC2=C3C(=C1)C4=CC=CC5=C4C(=CC=C5)C3=CC=C2"),
    ("picene", "C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC=CC=C54"),
    ("pleiadene", "C1=CC2=CC3=CC=CC4=C3C(=CC=C4)C=C2C=C1"),
    ("chrysene", "C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43"),
    ("pyrene", "C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2"),
    ("fluoranthene", "C1=CC=C2C(=C1)C3=CC=CC4=C3C2=CC=C4"),
    ("anthracene", "C1=CC=C2C=C3C=CC=CC3=CC2=C1"),
    ("phenanthrene", "C1=CC=C2C(=C1)C=CC3=CC=CC=C32"),
    ("phenalene", "C1C=CC2=CC=CC3=C2C1=CC=C3"),
    ("fluorene", "C1C2=CC=CC=C2C3=CC=CC=C31"),
    ("s-indacene", "C1=CC2=CC3=CC=CC3=CC2=C1"),
    ("as-indacene", "C1=CC2=C3C=CC=C3C=CC2=C1"),
    ("azulene", "C1=CC=C2C=CC=C2C=C1"),
    ("naphthalene", "C1=CC=C2C=CC=CC2=C1"),
    ("indene", "C1C=CC2=CC=CC=C21"),
    ("biphenylene", "C1=CC=C2C(=C1)C3=CC=CC=C23"),
    ("acenaphthylene", "C1=CC2=C3C(=C1)C=CC3=CC=C2"),
    ("aceanthrylene", "C1=CC=C2C3=C4C(=CC=CC4=CC2=C1)C=C3"),
    ("acephenanthrylene", "C1=CC=C2C(=C1)C=C3C=CC4=C3C2=CC=C4"),
]

_HETEROCYCLES = [
    ("phenazine", "C1=CC=C2C(=C1)N=C3C=CC=CC3=N2"),
    ("perimidine", "C1=CC2=C3C(=C1)NC=NC3=CC=C2"),
    ("acridine", "C1=CC=C2C(=C1)C=C3C=CC=CC3=N2"),
    ("phenanthridine", "C1=CC=C2C(=C1)C=NC3=CC=CC=C23"),
    ("carbazole", "C1=CC=C2C(=C1)C3=CC=CC=C3N2"),
    ("pteridine", "C1=CN=C2C(=N1)C=NC=N2"),
    ("cinnoline", "C1=CC=C2C(=C1)C=CN=N2"),
    ("quinazoline", "C1=CC=C2C(=C1)C=NC=N2"),
    ("quinoxaline", "C1=CC=C2C(=C1)N=CC=N2"),
    ("phthalazine", "C1=CC=C2C=NN=CC2=C1"),
    ("quinoline", "C1=CC=C2C(=C1)C=CC=N2"),
    ("isoquinoline", "C1=CC=C2C=NC=CC2=C1"),
    ("quinolizine", "C1C=CC=C2N1C=CC=C2"),
    ("purine", "C1=C2C(=NC=N1)N=CN2"),
    ("indazole", "C1=CC=C2C(=C1)C=NN2"),
    ("indole", "C1=CC=C2C(=C1)C=CN2"),
    ("isoindole", "C1=CC2=CNC=C2C=C1"),
    ("indolizine", "C1=CC2=CC=CN2C=C1"),
    ("pyrrolizine", "C1C=CN2C1=CC=C2"),
    ("xanthene", "C1C2=CC=CC=C2OC3=CC=CC=C31"),
    ("thioxanthene", "C1C2=CC=CC=C2SC3=CC=CC=C31"),
    ("selenoxanthene", "C1C2=CC=CC=C2[Se]C3=CC=CC=C31"),
    ("telluroxanthene", "C1C2=CC=CC=C2[Te]C3=CC=CC=C31"),
]
_ANTHRENES = {"O": "oxanthrene", "S": "thianthrene", "Se": "selenanthrene", "Te": "telluranthrene", "P": "phosphanthrene",
              "As": "arsanthrene", "B": "boranthrene", "Si": "silanthrene"}
_PHENO = {
    ("O", "N"): "phenoxazine", ("S", "N"): "phenothiazine", ("Se", "N"): "phenoselenazine", ("Te", "N"): "phenotellurazine",
    ("P", "N"): "phenazaphosphinine", ("As", "N"): "phenazarsinine", ("O", "S"): "phenoxathiine", ("O", "Se"): "phenoxaselenine",
    ("O", "Te"): "phenoxatellurine", ("O", "P"): "phenoxaphosphinine", ("O", "As"): "phenoxarsinine", ("O", "Sb"): "phenoxastibinine",
    ("S", "As"): "phenothiarsinine",
}
_P_AS_ANALOGUES = {
    "acridine": ("acridophosphine", "acridarsine"), "indole": ("phosphindole", "arsindole"),
    "indolizine": ("phosphindolizine", "arsindolizine"), "isoindole": ("isophosphindole", "isoarsindole"),
    "isoquinoline": ("isophosphinoline", "isoarsinoline"), "phenanthridine": ("phosphanthridine", "arsanthridine"),
    "quinoline": ("phosphinoline", "arsinoline"), "quinolizine": ("phosphinolizine", "arsinolizine"),
}

# retained numberings that do not follow P-25.3.3 (P-25.3.3): locants along the periphery, hetero positions, interior bonds
EXCEPTIONS = {
    "anthracene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {}, (("4a", "9a"), ("10a", "8a"))),
    "phenanthrene": (("1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "10", "10a"), {}, (("4a", "10a"), ("4b", "8a"))),
    "acridine": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "N"}, (("4a", "9a"), ("10a", "8a"))),
    "carbazole": (("1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "9a"), {"9": "N"}, (("4a", "9a"), ("4b", "8a"))),
    "xanthene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "O"}, (("4a", "9a"), ("10a", "8a"))),
    "thioxanthene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "S"}, (("4a", "9a"), ("10a", "8a"))),
    "selenoxanthene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "Se"}, (("4a", "9a"), ("10a", "8a"))),
    "telluroxanthene": (("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a"), {"10": "Te"}, (("4a", "9a"), ("10a", "8a"))),
    "purine": (("1", "2", "3", "4", "5", "6", "7", "8", "9"), {"1": "N", "3": "N", "7": "N", "9": "N"}, None),
    "cyclopenta[a]phenanthrene": (
        ("1", "2", "3", "4", "5", "6", "7", "8", "14", "15", "16", "17", "13", "12", "11", "9", "10"),
        {},
        (("5", "10"), ("8", "9"), ("13", "14")),
    ),
}
_PURINE_EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (3, 8), (8, 7), (7, 6), (6, 4)]


def _automorphic_images(graph, base):
    images, seen = [], set()
    for perm in graph.GetSubstructMatches(graph, uniquify=False, useChirality=False):
        image = {perm[a]: loc for a, loc in base.items()}
        ident = tuple(sorted(image.items()))
        if ident not in seen:
            seen.add(ident)
            images.append(image)
    return images


def exception_numbering(graph, name):
    locants, hetero, extra = EXCEPTIONS[name]
    n = len(locants)
    if extra is None:
        edges = _PURINE_EDGES
    else:
        index = {loc: i for i, loc in enumerate(locants)}
        edges = [(i, (i + 1) % n) for i in range(n)] + [(index[a], index[b]) for a, b in extra]
    query = skeleton([hetero.get(loc, "C") for loc in locants], edges)
    if query.GetNumAtoms() != graph.GetNumAtoms() or query.GetNumBonds() != graph.GetNumBonds():
        return None
    match = graph.GetSubstructMatch(query)
    if not match:
        return None
    return _automorphic_images(graph, {match[i]: loc for i, loc in enumerate(locants)})


def system_numberings(graph, name=None):
    """All numberings of a fused graph; the retained numberings of P-25.3.3 apply to the named exceptions."""
    graph = skeleton_of(graph)
    if name in EXCEPTIONS:
        images = exception_numbering(graph, name)
        if images:
            return images
    return fused_numberings(graph)


def _prefix_of(name):
    return name[:-1] + "o" if name.endswith("e") else name + "o"


@lru_cache(maxsize=None)
def _catalog():
    table = {}

    def add(name, prefix, kind, graph, **kw):
        table.setdefault(graph_key(graph), []).append(Prototype(name, prefix, kind, graph, exception=name if name in EXCEPTIONS else None, **kw))

    for name, smiles in _HYDROCARBONS:
        prefix = {"naphthalene": "naphtho", "anthracene": "anthra", "phenanthrene": "phenanthro"}.get(name, _prefix_of(name))
        add(name, prefix, "hydro", _from_smiles(smiles))
    for n in range(4, 13):
        term = numerical_term(n)
        add(_elide(term, "acene"), _elide(term, "acene")[:-1] + "o", "hydro", chain_graph([6] * n, [3] * (n - 2)))
        long_arm = (n + 1) // 2 if n % 2 else n // 2 + 1
        offsets = [3] * (long_arm - 2) + [2] + [3] * (n - long_arm - 1)
        add(term + "phene", term + "pheno", "hydro", chain_graph([6] * n, offsets))
    for n in range(6, 13):
        add(numerical_term(n) + "helicene", numerical_term(n) + "heliceno", "hydro", chain_graph([6] * n, [2] * (n - 2)))
    for n in (5, 7, 8, 9, 10, 11, 12):
        name = _elide(numerical_term(n), "alene")
        add(name, name[:-1] + "o", "hydro", chain_graph([n, n], []))
    for k in range(3, 7):
        add(numerical_term(k) + "phenylene", numerical_term(k) + "phenyleno", "hydro", _central_ring_graph(2 * k, [[6]] * k))
        add(numerical_term(k) + "naphthylene", numerical_term(k) + "naphthyleno", "hydro", _central_ring_graph(2 * k, [[6, 6]] * k))
    for name, smiles in _HETEROCYCLES:
        graph = _from_smiles(smiles)
        add(name, _prefix_of(name), "hetero", graph)
        if name in _P_AS_ANALOGUES:
            for element, analogue in zip(("P", "As"), _P_AS_ANALOGUES[name]):
                swapped = skeleton(
                    [element if a.GetSymbol() == "N" else a.GetSymbol() for a in graph.GetAtoms()],
                    [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in graph.GetBonds()],
                )
                add(analogue, _prefix_of(analogue), "hetero", swapped)
    base = _from_smiles("C1=CC=C2C(=C1)[Li]C3=CC=CC=C3[Be]2")
    marks = [i for i, a in enumerate(base.GetAtoms()) if a.GetSymbol() in ("Li", "Be")]
    pairs = {(e, e): n for e, n in _ANTHRENES.items()}
    pairs.update(_PHENO)
    for (first, second), name in pairs.items():
        symbols = [a.GetSymbol() for a in base.GetAtoms()]
        symbols[marks[0]], symbols[marks[1]] = first, second
        add(name, _prefix_of(name), "hetero", skeleton(symbols, [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in base.GetBonds()]))
    return table


def _images(proto, sub):
    """Numberings of `sub` (atom indices of the subgraph) given by the prototype's graph, or None if different."""
    if proto.fixed_numbering is not None:
        base = proto.fixed_numbering
        matches = sub.GetSubstructMatches(proto.graph, uniquify=False, useChirality=False)
        return [{m[a]: loc for a, loc in base.items()} for m in matches]
    numberings = system_numberings(proto.graph, proto.exception)
    if not numberings:
        raise UnsupportedStructure(f"no numbering for the component {proto.name}")
    base = numberings[0]
    return _unique([{m[a]: loc for a, loc in base.items()} for m in sub.GetSubstructMatches(proto.graph, uniquify=False, useChirality=False)])


def _unique(numberings):
    seen, out = set(), []
    for n in numberings:
        ident = tuple(sorted(n.items()))
        if ident not in seen:
            seen.add(ident)
            out.append(n)
    return out


def _cycle_order(sub):
    start = 0
    order = [start]
    previous = None
    current = start
    while True:
        nxt = [n.GetIdx() for n in sub.GetAtomWithIdx(current).GetNeighbors() if n.GetIdx() != previous]
        following = nxt[0]
        if following == start:
            return order
        order.append(following)
        previous, current = current, following


def _hw_stem(size, last, counts, saturated=False):
    if size == 6:
        return "ine" if last in _SIX_A | _SIX_B else "inine"
    if size == 3 and set(counts) == {"N"}:
        return "irine"
    return _STEMS[size]


def _join_hw(pieces):
    text = ""
    for piece in pieces:
        if text and piece[0] in "aeiou" and text.endswith("a"):
            text = text[:-1]
        text += piece
    return text


def monocycle_prototype(sub):
    """Prototype for a mancude monocycle (benzene, [n]annulene, Hantzsch-Widman and 'ine' names), numbering inside."""
    order = _cycle_order(sub)
    elements = [sub.GetAtomWithIdx(i).GetSymbol() for i in order]
    size = len(order)
    hetero = [e for e in elements if e != "C"]
    candidates = []
    for direction in (1, -1):
        for start in range(size):
            seq = [order[(start + direction * k) % size] for k in range(size)]
            locants = {atom: k + 1 for k, atom in enumerate(seq)}
            hetero_positions = [(locants[a], sub.GetAtomWithIdx(a).GetSymbol()) for a in seq if sub.GetAtomWithIdx(a).GetSymbol() != "C"]
            by_seniority = sorted(hetero_positions, key=lambda t: (HETERO_RANK[t[1]], t[0]))
            key = (sorted(p for p, _ in hetero_positions), [p for p, _ in by_seniority])
            candidates.append((key, seq))
    best_key = min(c[0] for c in candidates)
    best = [c[1] for c in candidates if c[0] == best_key]
    numberings = _unique([{atom: str(k + 1) for k, atom in enumerate(seq)} for seq in best])
    if not hetero:
        if size == 6:
            return "benzene", "benzo", numberings, "", True
        if size >= 7:
            return f"[{size}]annulene", "cyclo" + numerical_term(size), numberings, "", True
        stems = {3: "cyclopropa", 4: "cyclobuta", 5: "cyclopenta"}
        return "cyclo" + stems[size][5:-1] + "ene", stems[size], numberings, "", False
    seq = best[0]
    positions = {sub.GetAtomWithIdx(a).GetSymbol(): [] for a in seq}
    for k, atom in enumerate(seq):
        symbol = sub.GetAtomWithIdx(atom).GetSymbol()
        if symbol != "C":
            positions[symbol].append(k + 1)
    positions.pop("C", None)
    counts = {e: len(v) for e, v in positions.items()}
    citation = sorted(counts, key=lambda e: HETERO_RANK[e])
    locant_text = ",".join(str(p) for e in citation for p in sorted(positions[e]))
    retained = _retained_monocycle(size, positions, hetero)
    if retained:
        return retained, _RETAINED_PREFIX[retained], numberings, "", True
    if size > 10:
        pieces = [("" if counts[e] == 1 else numerical_term(counts[e])) + A_PREFIX[e] for e in citation]
        stem = "cyclo" + _alkane_stem(size) + "ine"
        name = "".join(pieces) + stem
        return name, _prefix_of(name), numberings, locant_text, True
    pieces = []
    for e in citation:
        if counts[e] > 1:
            pieces.append(numerical_term(counts[e]))
        pieces.append(A_PREFIX[e])
    stem = _hw_stem(size, citation[-1], counts)
    name = _join_hw(pieces + [stem])
    text = locant_text if len(hetero) > 1 else ""
    return name, _prefix_of(name), numberings, text, True


def _alkane_stem(size):
    from ._numerals import alkane_name

    return alkane_name(size)[:-3]


def _retained_monocycle(size, positions, hetero):
    key = (size, tuple(sorted(hetero)))
    if key == (5, ("N", "N")):
        return "pyrazole" if sorted(positions["N"]) == [1, 2] else "imidazole" if sorted(positions["N"]) == [1, 3] else None
    if key == (6, ("N", "N")):
        return {(1, 2): "pyridazine", (1, 3): "pyrimidine", (1, 4): "pyrazine"}.get(tuple(sorted(positions["N"])))
    return _RETAINED_MONOCYCLES.get(key)


def benzo_unit_prototype(sub, rings):
    """A benzene ring fused to a heteromonocycle of five or more members (P-25.2.2.4), treated as one component."""
    sizes = sorted(len(r) for r in rings)
    if len(rings) != 2 or sizes[0] != 6 and 6 not in sizes:
        return None
    hetero_ring = next((r for r in rings if any(sub.GetAtomWithIdx(a).GetSymbol() != "C" for a in r)), None)
    benzene = next((r for r in rings if all(sub.GetAtomWithIdx(a).GetSymbol() == "C" for a in r) and len(r) == 6), None)
    if hetero_ring is None or benzene is None or len(hetero_ring) < 5 or hetero_ring is benzene:
        return None
    shared = set(hetero_ring) & set(benzene)
    if len(shared) != 2:
        return None
    ring_sub = Chem.RWMol()
    mapping = {}
    for a in hetero_ring:
        mapping[a] = ring_sub.AddAtom(Chem.Atom(sub.GetAtomWithIdx(a).GetSymbol()))
    for k in range(len(hetero_ring)):
        ring_sub.AddBond(mapping[hetero_ring[k]], mapping[hetero_ring[(k + 1) % len(hetero_ring)]], Chem.BondType.SINGLE)
    mono = skeleton_of(ring_sub.GetMol())
    name, prefix, _, text, _ = monocycle_prototype(mono)
    candidates = []
    for ring_dir in (1, -1):
        order = list(hetero_ring)
        if ring_dir == -1:
            order.reverse()
        for start in range(len(order)):
            seq = order[start:] + order[:start]
            # locant 1: atom of the hetero ring next to a fusion atom, numbering runs away from the fusion bond
            if seq[0] in shared or seq[-1] not in shared:
                continue
            nonfusion = [a for a in seq if a not in shared]
            locants = {a: k + 1 for k, a in enumerate(nonfusion)}
            positions = [(locants[a], sub.GetAtomWithIdx(a).GetSymbol()) for a in nonfusion if sub.GetAtomWithIdx(a).GetSymbol() != "C"]
            by_seniority = sorted(positions, key=lambda t: (HETERO_RANK[t[1]], t[0]))
            candidates.append(((sorted(p for p, _ in positions), [p for p, _ in by_seniority]), seq, locants))
    if not candidates:
        return None
    best = min(c[0] for c in candidates)
    chosen = [c for c in candidates if c[0] == best][0]
    seq, locants = chosen[1], chosen[2]
    counts = {}
    for k, a in enumerate(seq):
        symbol = sub.GetAtomWithIdx(a).GetSymbol()
        if a not in shared and symbol != "C":
            counts.setdefault(symbol, []).append(locants[a])
    citation = sorted(counts, key=lambda e: HETERO_RANK[e])
    locant_text = ",".join(str(p) for e in citation for p in sorted(counts[e]))
    return name, prefix, locant_text


def identify(sub, rings=None):
    """Components matching the graph `sub` (a connected ring subset), as `Component` objects over its atom indices."""
    n = sub.GetNumAtoms()
    key = graph_key(sub)
    found = []
    for proto in _catalog().get(key, []):
        numberings = _images(proto, sub)
        found.append(_bind(proto, sub, numberings))
    ring_atoms = [list(r) for r in sub.GetRingInfo().AtomRings()]
    if len(ring_atoms) == 1:
        name, prefix, numberings, text, parent_ok = monocycle_prototype(sub)
        comp = Component(name, prefix, "mono" if parent_ok else "attached_only", list(range(n)), numberings,
                         [a.GetSymbol() for a in sub.GetAtoms()], [n], hetero_text=text)
        found.append(_finish(comp, sub))
    elif len(ring_atoms) == 2 and not found:
        unit = _benzo_unit(sub, ring_atoms)
        if unit is not None:
            found.append(unit)
    families = _family_components(sub, ring_atoms)
    found.extend(families)
    return found


def _bind(proto, sub, numberings):
    comp = Component(
        proto.name, proto.prefix, proto.kind, list(range(sub.GetNumAtoms())), numberings,
        [a.GetSymbol() for a in sub.GetAtoms()], sorted((len(r) for r in sub.GetRingInfo().AtomRings()), reverse=True),
        retained=proto.retained, benzo_unit=proto.benzo_unit, hetero_text=proto.hetero_text,
    )
    return _finish(comp, sub)


def _finish(comp, sub):
    comp.ring_sizes = sorted((len(r) for r in sub.GetRingInfo().AtomRings()), reverse=True)
    if comp.kind in ("mono", "attached_only"):
        comp.orientation = (-1, 0, 0, 0)
    else:
        try:
            comp.orientation = orientation_key(sub)
        except UnsupportedStructure:
            comp.orientation = (0, 0, 0, 0)
    comp.senior_key = seniority_key(comp, sub)
    return comp


def seniority_key(comp, sub):
    """P-25.3.2.4 (a)-(j): smaller is more senior."""
    hetero = [e for e in comp.elements if e != "C"]
    counts = {e: hetero.count(e) for e in set(hetero)}
    numbering = comp.numberings[0]
    by_atom = {a: numbering[a] for a in numbering}
    hetero_atoms = [i for i, e in enumerate(comp.elements) if e != "C"]
    hetero_locants = sorted(_locant_key(by_atom[i]) for i in hetero_atoms)
    ordered_locants = []
    for e in HETERO_ORDER:
        ordered_locants.extend(sorted(_locant_key(by_atom[i]) for i in hetero_atoms if comp.elements[i] == e))
    fusion_carbons = sorted(
        _locant_key(by_atom[i]) for i, e in enumerate(comp.elements) if e == "C" and sub.GetAtomWithIdx(i).GetDegree() >= 3 and by_atom[i][-1].isalpha()
    )
    return (
        min((PARENT_HETERO_RANK[e] for e in hetero), default=99),
        -len(comp.ring_sizes),
        tuple(-s for s in comp.ring_sizes),
        -len(hetero),
        -len(counts),
        tuple(-counts.get(e, 0) for e in HETERO_ORDER),
        comp.orientation,
        hetero_locants,
        ordered_locants,
        fusion_carbons,
    )


def _benzo_unit(sub, rings):
    unit = benzo_unit_prototype(sub, rings)
    if unit is None:
        return None
    name, hetero_name_prefix, locant_text = unit
    hetero_ring = next(r for r in rings if any(sub.GetAtomWithIdx(a).GetSymbol() != "C" for a in r))
    base_name = name
    glued = ("benz" if base_name[0] in "aeiou" else "benzo") + base_name
    if hetero_name_prefix == "imidazo":
        hetero_name_prefix = "imidazolo"
    prefix = ("benz" if hetero_name_prefix[0] in "aeiou" else "benzo") + hetero_name_prefix
    shared_numberings = _benzo_numberings(sub, rings, hetero_ring)
    comp = Component(glued, prefix, "hetero", list(range(sub.GetNumAtoms())), shared_numberings,
                     [a.GetSymbol() for a in sub.GetAtoms()], [], retained=False, benzo_unit=True, hetero_text=locant_text)
    return _finish(comp, sub)


def _benzo_numberings(sub, rings, hetero_ring):
    numberings = fused_numberings(sub)
    hetero_atoms = [a.GetIdx() for a in sub.GetAtoms() if a.GetSymbol() != "C"]
    keep = []
    for n in numberings:
        hetero_set = sorted(_locant_key(n[a]) for a in hetero_atoms)
        keep.append((hetero_set, n))
    low = min(k for k, _ in keep)
    return [n for k, n in keep if k == low]


def _family_components(sub, ring_atoms):
    """Phenanthroline and naphthyridine isomers (Table 2.8): the all-carbon skeleton with one nitrogen in each end ring."""
    out = []
    elements = [a.GetSymbol() for a in sub.GetAtoms()]
    if sorted(set(elements)) != ["C", "N"] or elements.count("N") != 2:
        return out
    plain = skeleton(["C"] * sub.GetNumAtoms(), [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in sub.GetBonds()])
    base = {graph_key(_from_smiles("C1=CC=C2C=CC=CC2=C1")): "naphthyridine", graph_key(_from_smiles("C1=CC=C2C(=C1)C=CC3=CC=CC=C32")): "phenanthroline"}
    name = base.get(graph_key(plain))
    if name is None:
        return out
    nitrogens = [i for i, e in enumerate(elements) if e == "N"]
    rings_of = {n: [k for k, r in enumerate(ring_atoms) if n in r] for n in nitrogens}
    if any(sub.GetAtomWithIdx(n).GetDegree() != 2 for n in nitrogens):
        return out
    if name == "phenanthroline":
        end_rings = [k for k, r in enumerate(ring_atoms) if sum(1 for other in ring_atoms if len(set(other) & set(r)) == 2) == 1]
        ok = len(end_rings) == 2 and all(any(k in end_rings for k in rings_of[n]) for n in nitrogens) and rings_of[nitrogens[0]] != rings_of[nitrogens[1]]
    else:
        ok = rings_of[nitrogens[0]] != rings_of[nitrogens[1]]
    if not ok:
        return out
    numberings = fused_numberings(sub)
    locs = sorted(_locant_key(numberings[0][n])[0] for n in nitrogens)
    comp = Component(name, _prefix_of(name), "hetero", list(range(sub.GetNumAtoms())), numberings, elements, [], hetero_text=f"{locs[0]},{locs[1]}")
    out.append(_finish(comp, sub))
    return out
