"""Axial descriptors of allenes and cumulenes read from the SMILES text (P-92.1.2.2.1, P-93.4.2.2, P-93.4.2.3).

RDKit drops extended tetrahedral (`[C@AL1]`) and cumulene double-bond stereo when it reads a SMILES string, so both are
read from the text of acyclic cumulenes. An even cumulene is an elongated tetrahedron whose four ligands are the
substituents of its two terminals: ranked by CIP priority (p1 > p2 near, d1 > d2 far), `@` on (p1, p2, d1, d2) is `P` and
`@@` is `M`. An odd cumulene is planar and takes E or Z from the higher-ranking substituent of each terminal. The
descriptor is cited with the lowest locant of the cumulated double bonds in the numbered chain of the finished name,
or at the front of the substituent that holds the axis (P-93.4.2.4).
"""

import re

from rdkit import Chem

from ._common import STEREO_OVERRIDE, InternalInvariantError, UnsupportedStructure
from ._numerals import numerical_term

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

    def close(self, opener, closer, kind, label):
        self.order[opener][self.order[opener].index(label)] = closer
        self.order[closer].append(opener)
        self.bond[(opener, closer)] = kind
        self.bond[(closer, opener)] = _FLIP.get(kind, kind)


def _parse(smiles, rings=False):
    """The atoms of a SMILES string with their neighbors in written order, a ring closure standing at its digit; None
    for a ring closure unless `rings`, and for a disconnected string."""
    graph, stack, previous, kind, i = _Graph(), [], None, "-", 0
    opened = {}
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
        elif (c.isdigit() or c == "%") and rings and previous is not None:
            label, i = (smiles[i : i + 3], i + 3) if c == "%" else (c, i + 1)
            if label in opened:
                opener, opener_kind = opened.pop(label)
                graph.close(opener, previous, kind if kind != "-" else opener_kind, label)
            else:
                opened[label] = (previous, kind)
                graph.order[previous].append(label)
            kind = "-"
        elif c.isdigit() or c in "%.":
            return None
        else:
            match = _ATOM.match(smiles, i)
            if match is None or match.group(0) == "[H]":
                return None
            previous = graph.add(match.group(0), previous, kind)
            kind = "-"
            i = match.end()
    return None if opened else graph


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
            found.append((first[:0:-1] + [centre] + second[1:], "P" if anticlockwise else "M"))
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
            found.append((chain, "Z" if sides[0] == sides[1] else "E"))
    return found


def _merge(name, token):
    existing = re.match(r"\(([^()]*)\)-", name)
    if existing is None:
        return f"({token})-{name}"
    tokens = existing.group(1).split(",")
    if any(not re.match(r"\d+[A-Za-z]?$", t) for t in tokens):
        return f"({token})-{name}"
    ordered = ",".join(sorted(tokens + [token], key=lambda t: int(re.match(r"\d+", t).group(0))))
    return f"({ordered})-{name[existing.end():]}"


_MULTIPLIER = r"(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca|bis|tris|tetrakis|pentakis|hexakis)"
_ENCLOSED = re.compile(r"[(\[{][^()\[\]{}]*[)\]}]")
_STEREO_PREFERENCE = (("Z",), ("R", "r", "M", "m"))
_PREFERRED = frozenset("ZRrMm")
_CARBONYL_ENDING = re.compile(r"(?:oate|oic|oyl|oxo|amide|imide|one|al)(?![a-z])")
_CITED_BEFORE = re.compile(r"\([^()]*[A-Za-z]\)-$")
_ENCLOSING_MARKS = ("()", "[]", "{}")


def _top_level(name):
    while True:
        stripped = _ENCLOSED.sub("", name)
        if stripped == name:
            return name
        name = stripped


def _chain_length(name):
    """Number of atoms of the chain that `name` is built on: the last numerical stem that stands before an
    ane/ene/yne ending, the longer stem winning where one ends in another ('undec', 'dec')."""
    top = _top_level(name)
    stems = {1: "meth", 2: "eth", 3: "prop", 4: "but"}
    best = None
    for n in range(1, 61):
        stem = stems.get(n) or numerical_term(n).removesuffix("a")
        for match in re.finditer(rf"{stem}(?:a-|an|en|yn|-)", top):
            key = (match.start() + len(stem), len(stem), n)
            best = key if best is None or key > best else best
    return best[2] if best else None


def _ene_locants(name):
    return sorted(int(n) for m in _ENE_LOCANTS.finditer(_top_level(name)) for n in m.group(1).split(","))


def _cited_locants(name, word):
    found = set()
    for match in re.finditer(re.escape(word), name):
        before = re.sub(rf"{_MULTIPLIER}$", "", re.sub(r"[(\[{]$", "", name[: match.start()]))
        listed = re.search(r"(?<![\w,])(\d+(?:,\d+)*)-$", before)
        if listed:
            found.update(int(n) for n in listed.group(1).split(","))
    return found


class _Scope:
    def __init__(self, mol):
        from ._common import adjacency, halogen_substituents

        self.mol = mol
        self.neighbors = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()}
        self.graph = adjacency(mol)
        self.halogens = halogen_substituents(mol)
        self.words = {}

    def word(self, root, attached_to):
        from ._substituents import name_branch

        if (root, attached_to) not in self.words:
            try:
                self.words[(root, attached_to)] = name_branch(self.graph, root, attached_to, self.halogens, mol=self.mol)[0]
            except UnsupportedStructure:
                self.words[(root, attached_to)] = None
        return self.words[(root, attached_to)]

    def side(self, root, attached_to):
        seen, stack = {attached_to}, [root]
        while stack:
            atom = stack.pop()
            if atom not in seen:
                seen.add(atom)
                stack.extend(self.neighbors[atom])
        return seen - {attached_to}

    def paths(self, region, through):
        for first in sorted(region):
            stack = [[first]]
            while stack:
                path = stack.pop()
                if len(path) > 1 and (through is None or through in path):
                    yield path
                stack.extend(path + [n] for n in self.neighbors[path[-1]] if n in region and n not in path)

    def bond_is_double(self, a, b):
        return self.mol.GetBondBetweenAtoms(a, b).GetBondType() == Chem.BondType.DOUBLE

    def carbonyl_agrees(self, path, name):
        on_path = set(path)
        oxo = any(
            self.mol.GetAtomWithIdx(other).GetAtomicNum() == 8 and self.bond_is_double(atom, other)
            for atom in path
            for other in self.neighbors[atom]
            if other not in on_path
        )
        return oxo == bool(_CARBONYL_ENDING.search(_top_level(name)))

    def fits(self, path, name):
        doubles = [i + 1 for i in range(len(path) - 1) if self.bond_is_double(path[i], path[i + 1])]
        if len(path) != _chain_length(name) or doubles != _ene_locants(name):
            return False
        on_path = set(path)
        for position, atom in enumerate(path, 1):
            for other in self.neighbors[atom]:
                word = None if other in on_path else self.word(other, atom)
                cited = _cited_locants(name, word) if word else ()
                if cited and position not in cited:
                    return False
        return True


def _leading_preference(word):
    group = re.match(r"\(([^()]*)\)-", word)
    return tuple(t[-1] not in _PREFERRED for t in group.group(1).split(",")) if group else ()


def _preference(tokens):
    return tuple(tuple(sorted(locant for locant, code in tokens if code in wanted)) for wanted in _STEREO_PREFERENCE)


def _depth(text):
    depth = deepest = 0
    for character in text:
        depth += (character in "([{") - (character in ")]}")
        deepest = max(deepest, depth)
    return deepest


def _enclosed(cited):
    opening, closing = _ENCLOSING_MARKS[_depth(cited) % 3]
    return opening + cited + closing


def _replace_branch(name, word, cited):
    """Puts `cited` where the substituent `word` stands, in enclosing marks one level above the deepest it holds
    (P-16.5.4)."""
    start = next((m.start() for m in re.finditer(re.escape(word), name) if not _CITED_BEFORE.search(name[: m.start()])), None)
    if start is None:
        return None
    end = start + len(word)
    enclosed = start > 0 and name[start - 1] in "([{" and name[end : end + 1] in ")]}"
    if not enclosed and (start == 0 or name[start - 1] == " ") and name[end : end + 1] in ("", " "):
        return name[:start] + cited + name[end:]
    return name[: start - enclosed] + _enclosed(cited) + name[end + enclosed :]


def _words(name):
    words, depth, start = [], 0, 0
    for index, character in enumerate(name):
        depth += (character in "([{") - (character in ")]}")
        if character == " " and depth == 0:
            words.append(name[start:index])
            start = index + 1
    return words + [name[start:]]


def _replace_branches(name, word, items):
    """Cites the substituent `word` once when every copy has the same descriptors, else each copy apart, in the order of
    the descriptors, where the multiplying prefix stood."""
    cited = {c for _, c, _ in items}
    if len(cited) == 1:
        return _replace_branch(name, word, next(iter(cited)))
    pattern = re.compile(rf"(?:(?P<locants>\d+(?:,\d+)*)-)?{_MULTIPLIER}[(\[{{]{re.escape(word)}[)\]}}]")
    match = pattern.search(name)
    if match is None:
        return None
    ordered = sorted(items, key=lambda item: (item[0] or 0, _leading_preference(item[1])))
    if match.group("locants") is None:
        pieces = " ".join(c for _, c, _ in sorted(items, key=lambda item: _leading_preference(item[1])))
    else:
        pieces = "-".join(f"{locant}-{_enclosed(c)}" for locant, c, _ in ordered)
    return name[: match.start()] + pieces + name[match.end() :]


def _cite_in(scope, axes, name, region, through):
    """(`name` with the descriptors of the `axes` it numbers, the axes left): an axis is cited at the locants of the chain
    that the name numbers (the path through `through` whose length, double bonds and substituent locants the name
    repeats), or inside the name of the substituent that holds it when the chain does not."""
    words = _words(name)
    if len(words) == 1:
        return _cite_word(scope, axes, name, region, through, False)
    remaining, original = list(axes), list(words)
    for index, word in enumerate(original):
        cited, left = _cite_word(scope, remaining, word, region, through, True)
        if len(left) < len(remaining):
            words[index], remaining = cited, left
    for index in reversed(range(len(original))):
        cited, left = _cite_word(scope, remaining, original[index], region, through, True)
        if len(left) < len(remaining):
            words[index : index + 1] = sorted([words[index], cited], key=_leading_preference)
            remaining = left
    return " ".join(words), remaining


def _cite_word(scope, axes, name, region, through, must_cover):
    fitting = [p for p in scope.paths(region, through) if scope.fits(p, name)]
    fitting = [p for p in fitting if scope.carbonyl_agrees(p, name)] or fitting
    covered = lambda path: [a for a in axes if set(a[0]) <= set(path)]
    path = ()
    if fitting and (not must_cover or any(covered(p) for p in fitting)):
        most = max(len(covered(p)) for p in fitting)
        options = {}
        for candidate in fitting:
            if len(covered(candidate)) == most:
                tokens = tuple(sorted((min(candidate.index(x) + 1 for x in run), code) for run, code in covered(candidate)))
                options.setdefault(tokens, candidate)
        tokens = min(options, key=_preference)
        path = options[tokens]
        for locant, code in tokens:
            name = _merge(name, f"{locant}{code}")
    remaining = [a for a in axes if not set(a[0]) <= set(path)]
    edges = [
        (atom, other)
        for atom in sorted(region)
        for other in scope.neighbors[atom]
        if other in region and other not in path
    ]
    placed = {}
    for atom, other in sorted(edges, key=lambda edge: -len(scope.side(edge[1], edge[0]))):
        branch = scope.side(other, atom)
        inside = [a for a in remaining if set(a[0]) <= branch] if branch <= region else []
        word = scope.word(other, atom) if inside else None
        if not word or word not in name or word == name:
            continue
        cited, left = _cite_in(scope, inside, word, branch, other)
        if not left:
            placed.setdefault(word, []).append((path.index(atom) + 1 if atom in path else None, cited, inside))
    for word, items in placed.items():
        replaced = _replace_branches(name, word, items)
        if replaced is not None:
            name = replaced
            consumed = [a for _, _, inside in items for a in inside]
            remaining = [a for a in remaining if a not in consumed]
    return name, remaining


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
    scope = _Scope(mol)
    atoms = {a.GetIdx() for a in mol.GetAtoms()}
    cited, left = _cite_in(scope, axes, name, atoms, None)
    if left:
        raise UnsupportedStructure("a cumulene axis whose atoms the name does not number")
    return cited


_AXIS_FROM_CENTRE = {"R": "M", "S": "P", "r": "m", "s": "p"}


def annotate_ring_allene_stereo(mol, smiles):
    """Writes the descriptors of an allene that has one ring atom as terminal on the atoms of `mol` (property
    `_StereoOverride`), since RDKit drops the allene mark: M/P on the ring terminal when its two ring ligands differ
    in rank, and when they tie (the ring is symmetric about the axis, so the ring centres and the axis are
    pseudoasymmetric together) the m/p of the axis and the r/s of the centres (P-93.5.1.4.2.2)."""
    if "@AL" not in smiles or not any(c.isdigit() for c in smiles):
        return
    graph = _parse(smiles, rings=True)
    if graph is None or len(graph.text) != mol.GetNumAtoms():
        return
    ranks = _ranks(smiles)
    for centre, text in enumerate(graph.text):
        tag = re.search(r"@AL([12])", text)
        if tag is None or len(graph.order[centre]) != 2:
            continue
        first, second = (_chain(graph, centre, side) for side in graph.order[centre])
        if len(first) != 2 or len(second) != 2:
            continue
        terminals = (first[-1], second[-1])
        in_ring = [mol.GetAtomWithIdx(t).IsInRing() for t in terminals]
        if sum(in_ring) != 1:
            continue
        try:
            written = [_ligands(graph, t, mol.GetAtomWithIdx(t).GetTotalNumHs()) for t in terminals]
        except UnsupportedStructure:
            continue
        ranked = [_by_priority(ligands, ranks) for ligands in written]
        ring_side = in_ring.index(True)
        if ranked[1 - ring_side] is None:
            continue
        if ranked[ring_side] is not None:
            swaps = sum(w != list(o) for w, o in zip(written, ranked))
            code = "P" if (tag.group(1) == "1") != bool(swaps % 2) else "M"
            mol.GetAtomWithIdx(terminals[ring_side]).SetProp(STEREO_OVERRIDE, code)
            continue
        labels = _pseudoasymmetric_labels(smiles, mol, centre, terminals, written, ring_side, tag.group(1))
        for atom_index, label in labels.items():
            mol.GetAtomWithIdx(atom_index).SetProp(STEREO_OVERRIDE, label)


def _pseudoasymmetric_labels(smiles, mol, centre, terminals, written, ring_side, tag):
    """{atom: label} of the molecule whose allene is replaced by the tetrahedral atom with the same four ligands and
    the same handedness (an allene is an elongated tetrahedron); the label of that atom is the descriptor of the
    axis. The two ring ligands of the ring terminal tie, so only the stereogenic units decide their order."""
    from ._common import assign_cip_labels
    from ._spiro_stereo import _parity

    raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if raw is None or raw.GetNumAtoms() != mol.GetNumAtoms():
        return {}
    ring, far = terminals[ring_side], terminals[1 - ring_side]
    editable = Chem.RWMol(raw)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    for ligand in written[1 - ring_side]:
        if ligand == "H":
            hydrogen = editable.AddAtom(Chem.Atom(1))
            editable.GetAtomWithIdx(hydrogen).SetIntProp("_orig", -1)
            editable.AddBond(ring, hydrogen, Chem.BondType.SINGLE)
        else:
            editable.AddBond(ring, ligand, Chem.BondType.SINGLE)
    for index in sorted((centre, far), reverse=True):
        editable.RemoveAtom(index)
    new_of = {a.GetIntProp("_orig"): a.GetIdx() for a in editable.GetAtoms()}
    hydrogen = new_of.get(-1)
    desired = [hydrogen if x == "H" else new_of[x] for x in written[0] + written[1]]
    merged = editable.GetAtomWithIdx(new_of[ring])
    merged.SetNoImplicit(True)
    merged.SetNumExplicitHs(0)
    actual = [n.GetIdx() for n in merged.GetNeighbors()]
    if sorted(actual) != sorted(desired):
        return {}
    anticlockwise = (tag == "1") != bool(_parity(actual, desired))
    merged.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW if anticlockwise else Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    surrogate = editable.GetMol()
    try:
        Chem.SanitizeMol(surrogate)
        surrogate = Chem.RemoveHs(surrogate)
        Chem.AssignStereochemistry(surrogate, cleanIt=True, force=True)
        assign_cip_labels(surrogate)
    except InternalInvariantError:
        raise
    except (ValueError, RuntimeError):
        return {}
    labels = {}
    for atom in surrogate.GetAtoms():
        if atom.HasProp("_CIPCode"):
            origin = atom.GetIntProp("_orig")
            code = atom.GetProp("_CIPCode")
            labels[origin] = _AXIS_FROM_CENTRE.get(code, code) if origin == ring else code
    return labels
