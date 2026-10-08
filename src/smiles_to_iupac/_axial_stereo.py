"""Axial descriptors of allenes and cumulenes read from the SMILES text (P-92.1.2.2.1, P-93.4.2.2, P-93.4.2.3).

RDKit drops extended tetrahedral (`[C@AL1]`) and cumulene double-bond stereo when it reads a SMILES string, so both are
read from the text of acyclic cumulenes. An even cumulene is an elongated tetrahedron whose four ligands are the
substituents of its two terminals: ranked by CIP priority (p1 > p2 near, d1 > d2 far), `@` on (p1, p2, d1, d2) is `P` and
`@@` is `M`. An odd cumulene is planar and takes E or Z from the higher-ranking substituent of each terminal. The
descriptor is cited with the lowest locant of the cumulated double bonds in the finished name.
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure

_ATOM = re.compile(r"\[[^\]]+\]|Cl|Br|[BCNOPSFIbcnops]")
_ENE_LOCANTS = re.compile(r"-(\d+(?:,\d+)*)-(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?en(?=[a-z-]|$)")
_FLIP = {"/": "\\", "\\": "/"}


class _Graph:
    def __init__(self):
        self.text = []
        self.order = []
        self.bond = {}

    def add(self, text, parent, kind):
        index = len(self.text)
        self.text.append(text)
        self.order.append([])
        if parent is not None:
            self.order[index].append(parent)
            self.order[parent].append(index)
            self.bond[(parent, index)] = kind
            self.bond[(index, parent)] = _FLIP.get(kind, kind)
        return index


def _parse(smiles):
    """The atoms of an acyclic SMILES string with their neighbors in written order; None otherwise."""
    graph, stack, previous, kind, i = _Graph(), [], None, "-", 0
    while i < len(smiles):
        c = smiles[i]
        if c in "-=#$:/\\":
            kind = c
            i += 1
        elif c == "(":
            stack.append(previous)
            i += 1
        elif c == ")":
            previous = stack.pop()
            i += 1
        elif c.isdigit() or c in "%.":
            return None
        else:
            match = _ATOM.match(smiles, i)
            if match is None or match.group(0) == "[H]":
                return None
            previous = graph.add(match.group(0), previous, kind)
            kind = "-"
            i = match.end()
    return graph


def _double_neighbors(graph, atom, excluded):
    return [n for n in graph.order[atom] if graph.bond[(atom, n)] == "=" and n != excluded]


def _chain(graph, start, towards):
    """Atoms of the cumulated chain that begins start=towards, as [start, towards, ...]."""
    chain = [start, towards]
    while True:
        following = _double_neighbors(graph, chain[-1], chain[-2])
        if not following:
            return chain
        chain.append(following[0])


def _ligands(graph, terminal, hydrogens):
    """Written order of the terminal's ligands, the implicit hydrogen following the preceding atom."""
    ligands = [n for n in graph.order[terminal] if graph.bond[(terminal, n)] != "="]
    if hydrogens:
        after_parent = bool(graph.order[terminal]) and graph.order[terminal][0] < terminal and graph.bond[(terminal, graph.order[terminal][0])] != "="
        ligands.insert(1 if after_parent else 0, "H")
    if len(ligands) != 2:
        raise UnsupportedStructure("a cumulene terminal without two ligands")
    return ligands


def _ranks(smiles):
    """Legacy CIP ranks of the atoms, read from a molecule parsed under the legacy stereo perception that sets them."""
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(True)
    try:
        ranked = Chem.MolFromSmiles(smiles)
        Chem.AssignStereochemistry(ranked, cleanIt=True, force=True)
        return {a.GetIdx(): a.GetIntProp("_CIPRank") if a.HasProp("_CIPRank") else 0 for a in ranked.GetAtoms()}
    finally:
        Chem.SetUseLegacyStereoPerception(legacy)


def _by_priority(ligands, ranks):
    """(higher, lower) of two ligands, a hydrogen below any atom; None when they tie."""
    rank = lambda x: -1 if x == "H" else ranks[x]
    first, second = ligands
    if rank(first) == rank(second):
        return None
    return (first, second) if rank(first) > rank(second) else (second, first)


def _even_axes(graph, mol, ranks):
    found = []
    for centre, text in enumerate(graph.text):
        tag = re.search(r"@AL([12])", text)
        if not tag:
            continue
        sides = graph.order[centre]
        if len(sides) != 2:
            raise UnsupportedStructure("an allene centre without two double bonds")
        first, second = _chain(graph, centre, sides[0]), _chain(graph, centre, sides[1])
        terminals = (first[-1], second[-1])
        written, ordered = [], []
        for terminal in terminals:
            ligands = _ligands(graph, terminal, mol.GetAtomWithIdx(terminal).GetTotalNumHs())
            ranked = _by_priority(ligands, ranks)
            if ranked is None:
                break
            written.append(ligands)
            ordered.append(ranked)
        else:
            swaps = sum(w != list(o) for w, o in zip(written, ordered))
            anticlockwise = (tag.group(1) == "1") != bool(swaps % 2)
            found.append((len(first) + len(second) - 2, "P" if anticlockwise else "M"))
    return found


def _odd_axes(graph, mol, ranks):
    found = []
    for start in range(len(graph.text)):
        if len(_double_neighbors(graph, start, None)) != 1:
            continue
        chain = _chain(graph, start, _double_neighbors(graph, start, None)[0])
        end = chain[-1]
        doubles = len(chain) - 1
        if end < start or doubles < 3 or doubles % 2 == 0 or len(_double_neighbors(graph, end, None)) != 1:
            continue
        sides = []
        for terminal in (start, end):
            ligands = _ligands(graph, terminal, mol.GetAtomWithIdx(terminal).GetTotalNumHs())
            ranked = _by_priority(ligands, ranks)
            marked = [n for n in ligands if n != "H" and graph.bond[(terminal, n)] in "/\\"]
            if ranked is None or not marked:
                sides = None
                break
            up = graph.bond[(terminal, marked[0])] == "/"
            sides.append(up if ranked[0] == marked[0] else not up)
        if sides:
            found.append((doubles, "Z" if sides[0] == sides[1] else "E"))
    return found


def _cumulated_start(name, double_bonds):
    runs = []
    for match in _ENE_LOCANTS.finditer(name):
        run = []
        for number in (int(x) for x in match.group(1).split(",")):
            if run and number != run[-1] + 1:
                runs.append(run)
                run = []
            run.append(number)
        runs.append(run)
    runs = [r for r in runs if len(r) == double_bonds]
    if len(runs) != 1:
        raise UnsupportedStructure("the cumulated double bonds are not cited with consecutive locants in the name")
    return runs[0][0]


def _merge(name, token):
    existing = re.match(r"\(([^()]*)\)-", name)
    if existing is None:
        return f"({token})-{name}"
    tokens = existing.group(1).split(",")
    if any(not re.match(r"\d+[A-Za-z]?$", t) for t in tokens):
        return f"({token})-{name}"
    ordered = ",".join(sorted(tokens + [token], key=lambda t: int(re.match(r"\d+", t).group(0))))
    return f"({ordered})-{name[existing.end():]}"


def cite_axial_stereo(smiles, name):
    if "@AL" not in smiles and not any(c in smiles for c in "/\\"):
        return name
    graph = _parse(smiles)
    mol = Chem.MolFromSmiles(smiles)
    if graph is None or mol is None or len(graph.text) != mol.GetNumAtoms():
        return name
    if not ("@AL" in smiles or "=C=C=" in smiles.replace("/", "").replace("\\", "")):
        return name
    ranks = _ranks(smiles)
    axes = _even_axes(graph, mol, ranks) + _odd_axes(graph, mol, ranks)
    if not axes:
        return name
    if len(axes) > 1:
        raise UnsupportedStructure("more than one cumulene axis")
    double_bonds, descriptor = axes[0]
    return _merge(name, f"{_cumulated_start(name, double_bonds)}{descriptor}")
