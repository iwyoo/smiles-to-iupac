"""Acyclic chain parents carrying any mix of supported groups (P-41, P-44.1.1,
P-44.3, P-45): the most senior class present becomes the suffix and every
other group, ring or branch is cited as a substituent prefix through
`name_branch` (P-29.3.3, P-29.4). Rings may only be substituents; a principal
group on a ring is not handled here.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    lowest_locant_set,
    multiplied_word,
    name_from_substituents,
    specified_stereo_elements,
    substituent_locant_set_and_citation,
)
from ._hetero_prefixes import is_functional_carbon
from ._multiplicative import _bare_key
from ._multiplicative_ring import _SUFFIX_WORDS, _citation_key, _join, _prefix_text, _suffix_text, monocycle_spec, numberings
from ._substituents import format_substituent_prefixes, name_branch

_SENIORITY = ["acid", "amide", "nitrile", "aldehyde", "ketone", "alcohol", "thiol", "amine"]
_TERMINAL = {"acid", "amide", "nitrile", "aldehyde"}
_MAX_ATOMS = 80


def _double_oxygens(mol, carbon):
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(carbon).GetNeighbors()
        if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]


def _single_neighbors(mol, carbon, atomic_num):
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(carbon).GetNeighbors()
        if n.GetAtomicNum() == atomic_num and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0
    ]


def _terminal_heteroatom(mol, idx, hydrogens):
    atom = mol.GetAtomWithIdx(idx)
    return atom.GetDegree() == 1 and atom.GetTotalNumHs() == hydrogens and not atom.GetFormalCharge()


def _group_of(mol, carbon):
    """(class, atoms owned by the group) for a principal-capable group on
    `carbon`, else None. Raises on carbon-bound groups this engine cannot
    name (esters, acid halides, ...)."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return None
    nitrogens = [
        n.GetIdx()
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 3.0
    ]
    if nitrogens:
        others = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogens[0]]
        if len(others) != 1 or others[0].GetAtomicNum() != 6:
            raise UnsupportedStructure("a cyanide not bonded to carbon is not a nitrile")
        return "nitrile", {nitrogens[0]}
    oxygens = _double_oxygens(mol, carbon)
    if oxygens:
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        hetero = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 6 and n.GetIdx() != oxygens[0]]
        if any(
            mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() != 1.0
            for n in atom.GetNeighbors()
            if n.GetIdx() != oxygens[0]
        ):
            raise UnsupportedStructure("a cumulated carbonyl (ketene) is not supported")
        if not hetero:
            if len(carbon_neighbors) == 1 and atom.GetTotalNumHs() == 1:
                return "aldehyde", {oxygens[0]}
            if len(carbon_neighbors) == 2:
                return "ketone", {oxygens[0]}
            raise UnsupportedStructure("a formaldehyde-type carbonyl is not named by the chain engine")
        if len(hetero) == 1:
            other = hetero[0]
            if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
                return "acid", {oxygens[0], other.GetIdx()}
            if other.GetAtomicNum() == 7 and _terminal_heteroatom(mol, other.GetIdx(), 2):
                return "amide", {oxygens[0], other.GetIdx()}
        return None
    for z, hydrogens, name in ((8, 1, "alcohol"), (16, 1, "thiol"), (7, 2, "amine")):
        for n in _single_neighbors(mol, carbon, z):
            if _terminal_heteroatom(mol, n, hydrogens):
                return name, {n}
    return None


def _paths(adj, eligible):
    ends = [a for a in eligible if sum(1 for n in adj[a] if n in eligible) <= 1] or list(eligible)
    paths = []
    for start in ends:
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            onward = [n for n in adj[node] if n in eligible and n not in path]
            if not onward:
                paths.append(path)
            for n in onward:
                stack.append((n, path + [n]))
    return paths


def name_polyfunctional(mol) -> str:
    _check_scope(mol)
    multiplicative = _multiplicative_name(mol)
    if multiplicative is not None:
        return multiplicative
    return _select(mol)[1]


def _check_scope(mol):
    if mol.GetNumAtoms() > _MAX_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("this molecule is out of scope for the polyfunctional chain engine")
    if specified_stereo_elements(mol):
        raise UnsupportedStructure("stereodescriptors are not supported by the polyfunctional chain engine yet")


def _select(mol, attach=None, n_names=()):
    """(key, name, parts) of the best parent. `attach`: an atom that carries a
    free valence (a multiplicative unit); it takes the lowest locant after
    the principal groups and multiple bonds."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())

    groups = {}
    for atom in mol.GetAtoms():
        found = _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], {})[atom.GetIdx()] = found[1]
    ring_groups = _ring_occurrences(mol)
    classes = set(groups) | {c for c, _, _ in ring_groups}
    principal = next((name for name in _SENIORITY if name in classes), None)
    for atom in mol.GetAtoms():
        if atom.GetIsotope() or atom.GetNumRadicalElectrons():
            raise UnsupportedStructure("isotopes and radicals are not supported by the polyfunctional chain engine")
        if atom.GetFormalCharge() and not _is_nitro_part(atom):
            raise UnsupportedStructure("charged atoms are not supported by the polyfunctional chain engine")
        if (
            atom.GetAtomicNum() == 6
            and not atom.IsInRing()
            and principal != "acid"
            and _is_ester_like(mol, atom.GetIdx())
        ):
            raise UnsupportedStructure("an ester outranks every parent this engine can build except an acid")
    if principal in (None, "amine") and any(_substituted_amine_nitrogen(mol, a) for a in mol.GetAtoms()):
        if attach is not None or n_names:
            raise UnsupportedStructure("N-substituted amines inside a unit are not handled by the chain engine")
        return _substituted_amine(mol, graph, halogens, aromatic_atoms, groups, ring_groups)

    if principal is None:
        if attach is not None:
            raise UnsupportedStructure("a unit without a principal group is not supported")
        return _plain_parent(mol, graph, halogens, aromatic_atoms)

    chain_best = None
    principal_atoms = groups.get(principal, {})
    if principal_atoms:
        owned = set().union(*principal_atoms.values())
        eligible = {
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6
            and not a.IsInRing()
            and (a.GetIdx() in principal_atoms or not is_functional_carbon(mol, a.GetIdx()))
        }
        for path in _paths(graph, eligible):
            for chain in (path, path[::-1]):
                candidate = _evaluate(
                    mol, graph, halogens, aromatic_atoms, chain, principal, principal_atoms, owned, attach, n_names
                )
                if chain_best is None or candidate[0] < chain_best[0]:
                    chain_best = candidate
        if chain_best is not None and chain_best[0][0] == 1:
            chain_best = None
    chain_count = -chain_best[0][0] if chain_best else 0

    ring_best = _best_ring(
        mol, graph, halogens, aromatic_atoms, principal, [g for g in ring_groups if g[0] == principal], n_names
    )
    ring_count = ring_best[0] if ring_best else 0

    if ring_count and ring_count >= chain_count:
        if attach is not None:
            raise UnsupportedStructure("a ring parent inside a multiplicative unit is not supported yet")
        return ring_best[1]
    if chain_best is None:
        raise UnsupportedStructure("no parent carries the principal group")
    if principal in _TERMINAL and chain_count < len(principal_atoms):
        raise UnsupportedStructure("principal groups that need a 'carbo' suffix are not supported yet")
    if principal in _TERMINAL and chain_best[0][1] == -1:
        raise UnsupportedStructure("one-carbon acid, amide, nitrile and aldehyde parents use retained names")
    return chain_best


def _plain_parent(mol, graph, halogens, aromatic_atoms):
    """Parent without a principal group: the monocycle when there is one
    (P-44.1.2.2), else the longest chain, with every substituent a prefix."""
    ring_info = mol.GetRingInfo()
    rings = [r for r in ring_info.AtomRings()]
    if rings:
        if len(rings) != 1 or any(ring_info.NumAtomRings(a) != 1 for a in rings[0]):
            raise UnsupportedStructure("several rings without a principal group are not named by the chain engine")
        spec = monocycle_spec(mol, rings[0])
        if spec is None:
            raise UnsupportedStructure("this ring is not supported as a parent by the chain engine")
        ring_set = set(rings[0])
        roots = [
            (r, n.GetIdx())
            for r in rings[0]
            for n in mol.GetAtomWithIdx(r).GetNeighbors()
            if n.GetIdx() not in ring_set
        ]
        if not roots:
            raise UnsupportedStructure("an unsubstituted ring is not a polyfunctional case")
        entries = [
            (r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots
        ]
        best = None
        for locants in numberings(spec):
            key = (
                tuple(sorted(locants[r] for r, _, _ in entries)),
                _citation_key([(locants[r], name) for r, name, _ in entries]),
            )
            if best is None or key < best[0]:
                best = (key, locants)
        if len(entries) == 1:
            _, only_name, only_compound = entries[0]
            prefix_text = format_substituent_prefixes(
                {only_name: {"locants": [1], "compound": only_compound}}, omit_locants=True
            )
        else:
            prefix_text = _prefix_text(entries, best[1])
        name = _join(prefix_text, spec.parent)
        return ((0,), name, None)

    eligible = {
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not is_functional_carbon(mol, a.GetIdx())
    }
    best = None
    for path in _paths(graph, eligible):
        for chain in (path, path[::-1]):
            candidate = _evaluate_plain(mol, graph, halogens, aromatic_atoms, chain)
            if best is None or candidate[0] < best[0]:
                best = candidate
    if best is None:
        raise UnsupportedStructure("no chain carbons")
    return best


def _evaluate_plain(mol, graph, halogens, aromatic_atoms, chain):
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    chain_set = set(chain)
    ene, yne = [], []
    for a, b in zip(chain, chain[1:]):
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order == 2.0:
            ene.append(position_of[a])
        elif order == 3.0:
            yne.append(position_of[a])
    entries = {}
    for atom in chain:
        for neighbor in graph[atom]:
            if neighbor in chain_set:
                continue
            name, compound = name_branch(graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            entries.setdefault(position_of[atom], []).append((name, compound))
    if not entries:
        raise UnsupportedStructure("an unsubstituted chain is not a polyfunctional case")
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
    length = len(chain)
    single = total_count == 1
    prefix = format_substituent_prefixes(grouped, omit_locants=length == 1 or (length == 2 and (ene or yne) and single))
    if length == 2 and (ene or yne):
        body = "ethene" if ene else "ethyne"
    else:
        body = name_from_substituents(length, ene, yne, "e")
    name = prefix + body
    key = (
        -length,
        -(len(ene) + len(yne)),
        -len(ene),
        lowest_locant_set(ene + yne),
        lowest_locant_set(ene),
        -total_count,
        locant_set,
        citation,
        name,
    )
    return key, name, None


_RING_SUFFIX = {
    "acid": "carboxylic_acid",
    "amide": "amide",
    "nitrile": "nitrile",
    "aldehyde": "aldehyde",
    "ketone": "ketone",
    "alcohol": "alcohol",
    "thiol": "thiol",
    "amine": "amine",
}


def _best_ring(mol, graph, halogens, aromatic_atoms, principal, occurrences, n_names=()):
    """(principal group count, (key, name, parts)) of the monocycle bearing
    the most principal groups, or None."""
    if not occurrences:
        return None
    ring_info = mol.GetRingInfo()
    candidates = []
    for ring in ring_info.AtomRings():
        here = [o for o in occurrences if o[1] in ring]
        if not here:
            continue
        if any(ring_info.NumAtomRings(a) != 1 for a in ring):
            raise UnsupportedStructure("a fused or bridged ring parent is not handled by the chain engine")
        spec = monocycle_spec(mol, ring)
        if spec is None:
            raise UnsupportedStructure("this ring is not supported as a parent by the chain engine")
        candidates.append((len(here), ring, spec, here))
    if not candidates:
        return None
    top = max(c[0] for c in candidates)
    leading = [c for c in candidates if c[0] == top]
    if len(leading) > 1:
        raise UnsupportedStructure("several rings bear the principal group; a multiplicative name is needed")
    count, ring, spec, here = leading[0]
    owned = set().union(*(o[2] for o in here))
    ring_set = set(ring)
    roots = [
        (r, n.GetIdx())
        for r in ring
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in ring_set and n.GetIdx() not in owned
    ]
    for _, n in roots:
        neighbor_ring = next((set(r) for r in ring_info.AtomRings() if n in r and not set(r) & ring_set), None)
        if neighbor_ring is not None and _bare_key(mol, neighbor_ring) == _bare_key(mol, ring_set):
            raise UnsupportedStructure("identical rings joined directly form a ring assembly (P-28)")
    entries = [
        (r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots
    ]
    principal_atoms = [o[1] for o in here]
    best = None
    for locants in numberings(spec):
        key = (
            tuple(sorted(locants[a] for a in principal_atoms)),
            tuple(sorted(locants[r] for r, _, _ in entries)),
            _citation_key([(locants[r], name) for r, name, _ in entries]),
        )
        if best is None or key < best[0]:
            best = (key, locants)
    locants = best[1]
    suffix_name = _RING_SUFFIX[principal]
    suffix_locants = [locants[a] for a in principal_atoms]
    prefix_text = _ring_prefix_text(entries, locants, n_names)
    if spec.kind == "cycloalkane" and count == 1 and not entries:
        word = _SUFFIX_WORDS[suffix_name]
        stem = spec.parent[:-1] if word[0] in "aeiouy" else spec.parent
        core = stem + word
    else:
        core, _ = _suffix_text(spec.parent, suffix_name, suffix_locants, spec)
    name = _join(prefix_text, core)
    return count, (((-count,), name, None), name, None)


def _substituted_amine(mol, graph, halogens, aromatic_atoms, groups, ring_groups):
    """A single secondary or tertiary amine: the ring or chain bound to
    nitrogen is the parent (ring first, P-44.1.2.2) and every other group on
    nitrogen is cited as an 'N-' prefix (P-62.2.2)."""
    amine_nitrogens = [a for a in mol.GetAtoms() if _substituted_amine_nitrogen(mol, a)]
    primaries = [a for a in groups.get("amine", {})] + [r for c, r, _ in ring_groups if c == "amine"]
    if len(amine_nitrogens) != 1 or primaries:
        raise UnsupportedStructure("several amine groups with N-substitution are not handled by the chain engine")
    nitrogen = amine_nitrogens[0]
    n_idx = nitrogen.GetIdx()
    if nitrogen.IsInRing() or nitrogen.GetTotalNumHs() > 1:
        raise UnsupportedStructure("a ring nitrogen is not an acyclic amine parent")
    neighbors = [n.GetIdx() for n in nitrogen.GetNeighbors()]
    if any(
        mol.GetAtomWithIdx(c).GetAtomicNum() != 6
        or _double_oxygens(mol, c)
        or is_functional_carbon(mol, c)
        or mol.GetBondBetweenAtoms(n_idx, c).GetBondTypeAsDouble() != 1.0
        for c in neighbors
    ):
        raise UnsupportedStructure("this nitrogen is not a plain amine nitrogen")
    arms = {c: _arm_atoms(graph, c, n_idx) for c in neighbors}
    if sum(len(a) for a in arms.values()) != len(set().union(*arms.values())) or n_idx in set().union(*arms.values()):
        raise UnsupportedStructure("a nitrogen closing a ring is not an acyclic amine parent")
    results = []
    for c in neighbors:
        others = [o for o in neighbors if o != c]
        n_names = [name_branch(graph, o, n_idx, halogens, aromatic_atoms, mol=mol, unsaturated=True) for o in others]
        parent, mapped = _amine_parent_molecule(mol, arms[c], c, n_idx)
        result = _select(parent, None, n_names)
        ring = parent.GetAtomWithIdx(mapped).IsInRing()
        results.append((ring, _ring_rank(parent, mapped), _chain_size(parent, mapped), result))
    results.sort(key=lambda r: (not r[0], tuple(-x for x in r[1]), -r[2], r[3][1]))
    return results[0][3]


def _amine_parent_molecule(mol, atoms, carbon, n_idx):
    editable = Chem.RWMol(mol)
    editable.GetAtomWithIdx(carbon).SetAtomMapNum(1)
    keep = set(atoms) | {n_idx}
    for idx in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(idx)
    parent = editable.GetMol()
    for a in parent.GetAtoms():
        if a.GetAtomicNum() == 7 and not a.IsInRing():
            a.SetNumExplicitHs(2)
            a.SetNoImplicit(True)
    Chem.SanitizeMol(parent)
    mapped = next(a.GetIdx() for a in parent.GetAtoms() if a.GetAtomMapNum() == 1)
    for a in parent.GetAtoms():
        a.SetAtomMapNum(0)
    return parent, mapped


def _ring_rank(mol, idx):
    ring = next((set(r) for r in mol.GetRingInfo().AtomRings() if idx in r), None)
    if ring is None:
        return (0, 0, 0)
    nitrogen = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 7 for a in ring)
    hetero = any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in ring)
    aromatic = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
    return (2 if nitrogen else 1 if hetero else 0, len(ring), int(aromatic))


def _chain_size(mol, idx):
    if mol.GetAtomWithIdx(idx).IsInRing():
        return 0
    graph = adjacency(mol)
    eligible = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
    return max((len(p) for p in _paths(graph, eligible) if idx in p), default=0)


_LINKER_WORDS = {8: "oxy", 16: "sulfanediyl", 34: "selanediyl"}


def _arm_atoms(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _multiplicative_name(mol):
    """P-15.3: identical chain parents joined through one oxygen, chalcogen,
    NH or N atom are named multiplicatively (2,2'-oxydi(ethan-1-ol))."""
    graph = adjacency(mol)
    for z in mol.GetAtoms():
        if z.IsInRing() or z.GetFormalCharge() or z.GetIsotope():
            continue
        neighbors = [n.GetIdx() for n in z.GetNeighbors()]
        hydrogens = z.GetTotalNumHs()
        number = z.GetAtomicNum()
        if number in _LINKER_WORDS and len(neighbors) == 2 and hydrogens == 0:
            linker, arms = _LINKER_WORDS[number], 2
        elif number == 7 and len(neighbors) == 2 and hydrogens == 1:
            linker, arms = "azanediyl", 2
        elif number == 7 and len(neighbors) == 3 and hydrogens == 0:
            linker, arms = "nitrilo", 3
        else:
            continue
        if any(
            mol.GetAtomWithIdx(n).GetAtomicNum() != 6
            or _double_oxygens(mol, n)
            or is_functional_carbon(mol, n)
            or mol.GetAtomWithIdx(n).GetIsAromatic()
            for n in neighbors
        ):
            continue
        parts_atoms = [_arm_atoms(graph, n, z.GetIdx()) for n in neighbors]
        if sum(len(p) for p in parts_atoms) != mol.GetNumAtoms() - 1 or any(
            parts_atoms[i] & parts_atoms[j] for i in range(arms) for j in range(i)
        ):
            continue
        units = [_unit_molecule(mol, atoms, n) for atoms, n in zip(parts_atoms, neighbors)]
        keys = {Chem.MolToSmiles(unit[0]) for unit in units}
        if len(keys) != 1:
            continue
        unit, attach = units[0]
        _, _, parts = _select(unit, attach)
        prefix, body, tail, locant = parts
        locants = ",".join(str(locant) + "'" * i for i in range(arms))
        text = prefix + body
        if prefix:
            word = {2: "bis", 3: "tris"}[arms]
            return f"{locants}-{linker}{word}({text}){tail}"
        word = {2: "di", 3: "tri"}[arms]
        unit_text = f"({text})" if any(ch.isdigit() or ch == "-" for ch in text) else text
        return f"{locants}-{linker}{word}{unit_text}{tail}"
    return None


def _unit_molecule(mol, atoms, attach):
    editable = Chem.RWMol(mol)
    editable.GetAtomWithIdx(attach).SetAtomMapNum(1)
    for idx in sorted(set(range(mol.GetNumAtoms())) - atoms, reverse=True):
        editable.RemoveAtom(idx)
    unit = editable.GetMol()
    Chem.SanitizeMol(unit)
    attach_idx = next(a.GetIdx() for a in unit.GetAtoms() if a.GetAtomMapNum() == 1)
    canonical = Chem.Mol(unit)
    for a in canonical.GetAtoms():
        a.SetAtomMapNum(1 if a.GetIdx() == attach_idx else 0)
    return canonical, attach_idx


def _is_nitro_part(atom):
    """The charged atoms of a nitro group: N+ bonded to two oxygens (one O-)."""
    if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        return len(oxygens) == 2 and sum(o.GetFormalCharge() for o in oxygens) == -1 and atom.GetDegree() == 3
    if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1 and atom.GetDegree() == 1:
        (n,) = atom.GetNeighbors()
        return _is_nitro_part(n)
    return False


def _is_ester_like(mol, carbon):
    atom = mol.GetAtomWithIdx(carbon)
    if not _double_oxygens(mol, carbon):
        return False
    return any(
        n.GetAtomicNum() in (8, 7, 16, 9, 17, 35, 53) and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0
        and not (n.GetAtomicNum() == 8 and _terminal_heteroatom(mol, n.GetIdx(), 1))
        and not (n.GetAtomicNum() == 7 and _terminal_heteroatom(mol, n.GetIdx(), 2))
        for n in atom.GetNeighbors()
    )


def _substituted_amine_nitrogen(mol, atom):
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() or atom.GetIsAromatic():
        return False
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return len(carbons) >= 2


def _ring_occurrences(mol):
    """[(class, ring_atom, owned atoms)] for every principal-capable group
    sitting directly on a ring atom (or, for a ketone, the ring carbonyl)."""
    found = []
    for atom in mol.GetAtoms():
        if not atom.IsInRing():
            continue
        r = atom.GetIdx()
        for n in atom.GetNeighbors():
            if n.IsInRing():
                continue
            z, i = n.GetAtomicNum(), n.GetIdx()
            order = mol.GetBondBetweenAtoms(r, i).GetBondTypeAsDouble()
            if z == 8 and order == 2.0 and n.GetDegree() == 1:
                found.append(("ketone", r, {i}))
            elif z == 8 and _terminal_heteroatom(mol, i, 1):
                found.append(("alcohol", r, {i}))
            elif z == 16 and _terminal_heteroatom(mol, i, 1):
                found.append(("thiol", r, {i}))
            elif z == 7 and _terminal_heteroatom(mol, i, 2):
                found.append(("amine", r, {i}))
            elif z == 6:
                group = _group_of(mol, i)
                if group is not None and group[0] in _TERMINAL:
                    found.append((group[0], r, {i} | group[1]))
    return found


def _with_n_names(grouped, n_names):
    if not n_names:
        return grouped
    merged = {name: {"locants": list(info["locants"]), "compound": info["compound"]} for name, info in grouped.items()}
    for name, compound in n_names:
        merged.setdefault(name, {"locants": [], "compound": compound})["locants"].append("N")
    return merged


def _ring_prefix_text(entries, locants, n_names):
    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locants[r])
    return format_substituent_prefixes(_with_n_names(grouped, n_names)) if (grouped or n_names) else ""


def _evaluate(mol, graph, halogens, aromatic_atoms, chain, principal, principal_atoms, owned, attach=None, n_names=()):
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    chain_set = set(chain)
    if attach is not None and attach not in chain_set:
        return ((1,), "", None)
    on_chain = [a for a in principal_atoms if a in chain_set]
    if not on_chain or (principal in _TERMINAL and any(position_of[a] not in (1, len(chain)) for a in on_chain)):
        return ((1,), "", None)
    ene, yne = [], []
    for a, b in zip(chain, chain[1:]):
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order == 2.0:
            ene.append(position_of[a])
        elif order == 3.0:
            yne.append(position_of[a])
    entries = {}
    for atom in chain:
        for neighbor in graph[atom]:
            if neighbor in chain_set or neighbor in owned:
                continue
            name, compound = name_branch(graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            entries.setdefault(position_of[atom], []).append((name, compound))
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
    suffix_locants = sorted(position_of[a] for a in on_chain)
    count = len(on_chain)
    length = len(chain)
    force = attach is not None
    prefix = format_substituent_prefixes(
        _with_n_names(grouped, n_names), omit_locants=length == 1 and not force and not n_names
    )
    tail = ""
    if principal == "acid":
        body, tail = name_from_substituents(length, ene, yne, multiplied_word(count, "oic")), " acid"
    elif principal == "amide":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "amide"))
    elif principal == "nitrile":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "nitrile"))
    elif principal == "aldehyde":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "al"))
    else:
        word = {"ketone": "one", "alcohol": "ol", "thiol": "thiol", "amine": "amine"}[principal]
        body = name_from_substituents(
            length, ene, yne, multiplied_word(count, word), suffix_locants, force_own_locant=force
        )
    name = prefix + body + tail
    attach_locant = position_of[attach] if attach is not None else 0
    key = (
        -count,
        -length,
        -(len(ene) + len(yne)),
        -len(ene),
        tuple(suffix_locants),
        lowest_locant_set(ene + yne),
        lowest_locant_set(ene),
        attach_locant,
        -total_count,
        locant_set,
        citation,
        name,
    )
    return key, name, (prefix, body, tail, attach_locant)
