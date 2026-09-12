"""Shared graph utilities and validation for saturated hydrocarbon parent
hydrides, used by both the acyclic (`_acyclic.py`) and monocyclic (`_cyclic.py`)
naming modules.

- P-35.2.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'fluoro',
  'chloro', 'bromo', and 'iodo' are the preselected substituent prefixes for
  -F, -Cl, -Br, and -I respectively. These halogen atoms are always
  monovalent, non-skeletal substituents (P-44.3's "skeletal atoms" always
  means carbon in this repository) — never part of a parent hydride's
  counted chain/ring — so `carbon_adjacency` below lets chain/ring-skeleton
  search ignore them while substituent detection still finds them.
- P-92 (Chapter P-9, https://iupac.qmul.ac.uk/BlueBook/P9.html):
  `specified_stereocenters` below delegates all CIP priority-rule
  computation (atomic number, duplicate-atom treatment, mass number,
  pseudoasymmetry, ...) to RDKit's
  `rdCIPLabeler` rather than reimplementing P-92's rules directly -- this
  project's own contribution is only formatting the resulting label(s)
  into a name, not computing them. P-91.3: when two or more are cited
  together, they're joined in ascending locant order, comma-separated,
  in one parenthesized group, e.g. '(2S,3S)-3-chloro-2-hydroxybutanoic
  acid' (a Blue Book worked example) -- confirmed directly from the
  primary source text, not derived from this project's own reasoning.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._numerals import numerical_term


class UnsupportedStructure(NotImplementedError):
    pass


HALOGEN_PREFIXES = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_ALLOWED_ATOMIC_NUMS = {6, *HALOGEN_PREFIXES}

ENE_BOND_ORDER = 2.0
YNE_BOND_ORDER = 3.0


def validate_atoms_and_bonds(mol):
    """Structure-independent checks shared by every parent hydride kind: a
    single-fragment, all-carbon skeleton optionally bearing monovalent
    halogen substituents (P-35.2.1). Bond order (all single, or exactly one
    double/triple bond) and ring shape (none, one simple ring, or more) are
    checked separately by each naming module, since what's allowed there
    differs (see `_acyclic.py`, `_cyclic.py`, `_unsaturated.py`)."""
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than halogen substituents (F, Cl, Br, I; "
                "see P-35.2.1) are not supported yet (see P-21.2.3, skeletal "
                "replacement nomenclature)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
        elif atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1); "
                "polyvalent or bridging halogen structures are not supported"
            )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute (see P-44.3)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )


def non_single_bonds(mol):
    """List of (begin_atom_idx, end_atom_idx, bond_order) for every bond whose
    order isn't 1.0 (single). Used to classify a molecule's degree of
    unsaturation for dispatch (see `core.py`)."""
    return [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0
    ]


def adjacency(mol):
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        graph[a].append(b)
        graph[b].append(a)
    return graph


def carbon_adjacency(mol):
    """Like `adjacency`, but restricted to carbon atoms and the bonds directly
    between them. Chain/ring-skeleton search (P-44.3.2) must use this instead
    of `adjacency` so a terminal halogen substituent is never mistaken for a
    chain-extending skeletal atom; substituent-detection code should keep
    using the full `adjacency(mol)` so it can still find that halogen."""
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def ring_chain_attachment(graph, ring_atoms, excluded):
    """If `ring_atoms` (a plain, otherwise-unsubstituted monocyclic ring)
    has exactly one exocyclic branch, hanging off exactly one ring atom,
    return (ring_atom, chain_root); else None. `excluded` are atom indices
    ignored as branch roots (e.g. a ring hydroxyl oxygen already accounted
    for separately) -- pass an empty set/frozenset if there are none."""
    attachments = []
    for atom in ring_atoms:
        branch_roots = [n for n in graph[atom] if n not in ring_atoms and n not in excluded]
        if branch_roots:
            attachments.append((atom, branch_roots))
    if len(attachments) != 1:
        return None
    ring_atom, branch_roots = attachments[0]
    if len(branch_roots) != 1:
        return None
    return ring_atom, branch_roots[0]


def ring_chain_attachment_with_halogens(graph, ring_atoms, excluded, halogens):
    """Like `ring_chain_attachment`, but tolerates the ring's other atoms
    each carrying a single plain halogen substituent (an exocyclic
    neighbor found in `halogens`) instead of requiring `ring_atoms` to be
    completely unsubstituted apart from the one chain -- lets a "phenyl
    chain" module recognize e.g. 4-chlorophenyl the same way it already
    recognizes plain phenyl (`name_branch`'s own halogenated-phenyl path
    then names the ring). Returns (ring_atom, chain_root) for the sole
    non-halogen exocyclic branch; else None (no such branch, more than one
    non-halogen branch, or any ring atom with more than one exocyclic
    neighbor)."""
    chain_attachment = None
    for atom in ring_atoms:
        branch_roots = [n for n in graph[atom] if n not in ring_atoms and n not in excluded]
        if not branch_roots:
            continue
        if len(branch_roots) != 1:
            return None
        (branch_root,) = branch_roots
        if branch_root in halogens:
            continue
        if chain_attachment is not None:
            return None
        chain_attachment = (atom, branch_root)
    return chain_attachment


def plain_saturated_ring_substituent_atoms(mol, graph, coming_from, root):
    """Ring atom set if `root` sits on a single plain, unsubstituted,
    saturated monocyclic ring whose only exocyclic bond is to
    `coming_from` -- else empty (a ring with a substituent/unsaturation, a
    polycyclic/spiro shape, or no ring at all). Used by any module naming
    a plain saturated ring as a substituent attached through one specific
    neighbor (e.g. `_carbamate.py`'s/`_ester.py`'s ring-on-the-alkoxy-side
    path) -- `name_branch`'s own `_simple_ring_substituent` already names
    the ring itself once it's let through the caller's own ring
    rejection; this only answers whether that shape applies here."""
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = set(ring)
        if root not in ring_atoms:
            continue
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
            continue
        attachment = ring_chain_attachment(graph, ring_atoms, set())
        if attachment is None:
            continue
        ring_atom, chain_root = attachment
        if ring_atom == root and chain_root == coming_from:
            return ring_atoms
    return set()


def is_plain_benzene_ring(mol, ring_atoms):
    """True if `ring_atoms` is exactly a 6-membered ring of aromatic carbons
    (a plain, unsubstituted-shape benzene ring) -- used by any chain-parent
    module (`_carboxylic_acid.py`, `_ketone.py`, ...) that names such a ring
    as a 'phenyl' substituent prefix on the chain rather than as the parent
    itself."""
    return len(ring_atoms) == 6 and all(
        mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic()
        for idx in ring_atoms
    )


_HETEROAROMATIC_MONOCYCLE_NAMES = {
    (6, 7): "pyridine",
    (5, 8): "furan",
    (5, 16): "thiophene",
    (5, 7): "pyrrole",
}


def heteroaromatic_monocycle_name(mol, ring_order):
    """Parent hydride name if `ring_order` (an ordered ring walk, e.g. from
    `ring_cycle`) is a plain, fully aromatic 5- or 6-membered monocycle
    with exactly one heteroatom matching one of the four simple
    heteroaromatic monocycles this project names as a "-yl" substituent
    prefix -- pyridine, furan, thiophene, or pyrrole (P-29.3.4.1) -- else
    None. Composition only (no charges/isotopes, exactly one heteroatom,
    the rest aromatic carbon); doesn't check for a single exocyclic
    attachment point itself -- see `ring_chain_attachment` for that."""
    n = len(ring_order)
    if n not in (5, 6):
        return None
    if any(
        mol.GetAtomWithIdx(idx).GetFormalCharge() != 0 or mol.GetAtomWithIdx(idx).GetIsotope() != 0
        for idx in ring_order
    ):
        return None
    if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_order):
        return None
    heteroatoms = [idx for idx in ring_order if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(heteroatoms) != 1:
        return None
    (heteroatom,) = heteroatoms
    return _HETEROAROMATIC_MONOCYCLE_NAMES.get((n, mol.GetAtomWithIdx(heteroatom).GetAtomicNum()))


def heteroaromatic_monocycle_yl_name(mol, ring_order, attachment_atom):
    """"-yl" substituent name (e.g. "pyridin-3-yl", "1H-pyrrol-2-yl") for a
    plain heteroaromatic monocycle recognized by
    `heteroaromatic_monocycle_name`, with the free valence at
    `attachment_atom` -- the heteroatom is fixed at locant 1 (its own
    established parent-hydride numbering), and `attachment_atom` gets
    whichever of the two ring-walk directions gives it the lower locant
    (P-29.3.4.1's own worked example, "pyridin-2-yl"). Pyrrole's own name
    always cites its indicated hydrogen ("1H-pyrrole" is itself the
    correct parent name, P-25.7.1.3) except when the substitution is
    directly at that N-H position, which consumes it instead (plain
    "pyrrol-1-yl", no citation needed, mirroring
    `_pyridine_heterocycle_fusion.py`'s identical `GetTotalNumHs() > 0`
    heuristic for the same tautomer distinction); pyridine/furan/thiophene
    never need this, since their heteroatom carries no H to begin with.
    Returns None if `ring_order` isn't one of the four recognized rings."""
    name = heteroaromatic_monocycle_name(mol, ring_order)
    if name is None:
        return None
    n = len(ring_order)
    (heteroatom,) = [idx for idx in ring_order if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    start = ring_order.index(heteroatom)
    target = ring_order.index(attachment_atom)
    forward = (target - start) % n
    backward = (start - target) % n
    locant = min(forward, backward) + 1
    stem = name[:-1] if name.endswith("e") else name
    indicated_hydrogen = ""
    if attachment_atom != heteroatom and mol.GetAtomWithIdx(heteroatom).GetTotalNumHs() > 0:
        indicated_hydrogen = "1H-"
    return f"{indicated_hydrogen}{stem}-{locant}-yl"


def plain_phenyl_substituent_atoms(mol, graph, roots):
    """Union of ring atoms for every plain, unsubstituted benzene ring in
    `mol` that hangs directly off one of `roots` with no other exocyclic
    attachment -- i.e. a lone 'phenyl' substituent directly on one of
    `roots`."""
    atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = set(ring)
        if not is_plain_benzene_ring(mol, ring_atoms):
            continue
        attachment = ring_chain_attachment(graph, ring_atoms, set())
        if attachment is None:
            continue
        ring_atom, _ = attachment
        if ring_atom in roots:
            atoms |= ring_atoms
    return atoms


def ordered_chain(graph, root, coming_from, excluded):
    """The chain of atoms starting at `root` and extending away from
    `coming_from`, ignoring `excluded` atoms (e.g. hydroxyl/carboxyl
    oxygens) the same way `_substituents._longest_chains_from_root`
    ignores halogens; None if it branches (more than one non-excluded,
    non-`previous` neighbor at any point)."""
    chain = [root]
    previous, current = coming_from, root
    while True:
        neighbors = [n for n in graph[current] if n != previous and n not in excluded]
        if not neighbors:
            return chain
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]
        chain.append(current)


def longest_branched_chain(graph, source, ring_boundary, excluded=frozenset(), halogens=frozenset()):
    """Generalizes `ordered_chain`: the longest simple chain starting at
    `source` (typically a suffix's own principal-characteristic-group
    carbon, always a tree leaf once its own heteroatoms are excluded)
    and extending through the tree, absorbing any branch encountered
    along the way into the chain itself whenever doing so makes it
    longer (P-44.3.2: the parent chain is one of the longest chains
    containing the principal characteristic group) -- e.g. ibuprofen's
    alpha-methyl becomes part of the parent chain ('propanoic acid'),
    not a separate '2-methyl-...' prefix on a shorter 'ethanoic acid'
    (PubChem CID 3672, '2-[4-(2-methylpropyl)phenyl]propanoic acid').

    `ring_boundary`: atoms the chain itself may never enter (typically
    every ring atom -- this project's own "phenyl chain" modules always
    cite a ring hanging off the chain as a substituent prefix, never as
    part of the parent chain), but which still shows up as an ordinary
    branch wherever it's adjacent to a chosen chain atom.
    `excluded`: atoms that are neither part of the chain nor ever cited
    as a separate branch (typically the suffix's own carbonyl/hydroxyl
    oxygens, already accounted for by the suffix name itself).
    `halogens`: halogen atom indices (typically `halogen_substituents(mol)`,
    the same dict passed to `name_branch`) -- excluded from the chain
    itself (a halogen is always monovalent and would otherwise dead-end
    the BFS at a leaf, silently absorbing it into the chain as if it were
    carbon), but, like `ring_boundary`, still surfaced as an ordinary
    branch wherever one is adjacent to a chosen chain atom (P-35.2.1: a
    halogen is always cited as a substituent prefix, never part of the
    parent chain).

    On a tie for longest, prefers the chain giving the greater number of
    substituents cited as prefixes (P-44.3.2's own next tie-break after
    chain length), then the lowest set of locants among those (P-14.5.2)
    -- e.g. 'c1ccccc1CC(C)C(=O)O' (both a phenyl-bearing carbon and a
    methyl-bearing carbon are one bond from the acid's C2, tied for
    farthest) is '2-methyl-3-phenylpropanoic acid' (two substituents:
    'methyl' + 'phenyl'), not '2-(phenylmethyl)propanoic acid' (one
    compound 'benzyl'-shaped substituent) -- PubChem PUG REST verified.

    Returns (chain, branches): `chain` is the winning path as an
    atom-index list with `source` first (so 1-based `enumerate(chain,
    start=1)` locants match the suffix's own fixed-C1 convention);
    `branches` is {position -> [branch_root_atom, ...]} for every
    neighbor of a chain atom that isn't itself on the chain or in
    `excluded` -- name each via `name_branch`, same as any other
    substituent (a ring-atom branch root is named as a ring substituent
    automatically, since `name_branch` already recognizes one)."""
    blocked = set(ring_boundary) | set(excluded) | set(halogens)
    dist = {source: 0}
    parent = {source: None}
    queue = [source]
    while queue:
        next_queue = []
        for node in queue:
            for neighbor in graph[node]:
                if neighbor in blocked or neighbor in dist:
                    continue
                dist[neighbor] = dist[node] + 1
                parent[neighbor] = node
                next_queue.append(neighbor)
        queue = next_queue
    farthest = max(dist.values())

    def branches_for(chain):
        chain_set = set(chain)
        branches = {}
        for position, atom in enumerate(chain, start=1):
            roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
            if roots:
                branches[position] = roots
        return branches

    best_chain, best_branches, best_key = None, None, None
    for node, d in dist.items():
        if d != farthest:
            continue
        chain = path_between(parent, source, node)
        branches = branches_for(chain)
        substituent_count = sum(len(roots) for roots in branches.values())
        locants = lowest_locant_set(pos for pos, roots in branches.items() for _ in roots)
        key = (-substituent_count, locants)
        if best_key is None or key < best_key:
            best_key, best_chain, best_branches = key, chain, branches
    return best_chain, best_branches


def longest_branched_chain_through(graph, required, ring_boundary, excluded=frozenset(), halogens=frozenset()):
    """Like `longest_branched_chain`, but `required` need not be a chain
    terminus (e.g. a ketone's own carbonyl carbon, always internal once
    its aryl-ketone case is separately rejected) -- finds one of the
    longest chains in the tree that includes `required` somewhere along
    it, by combining `required`'s two longest 'arms' (P-44.3.2, the same
    "longest chain containing the principal characteristic group" rule,
    here allowing the group to sit anywhere on the chain instead of
    fixing it at C1).

    On a tie for longest, prefers the combination giving the greater
    number of substituents cited as prefixes, then the lowest set of
    locants among those -- the same P-44.3.2 tie-break `longest_branched_
    chain` already applies, generalized to the two-arm case (found via
    real-data testing: a branch point with a halogen-bearing 1-carbon arm
    tied in length against a plain methyl arm, e.g. 'C(F)(F)Br' vs 'C'
    off the same alcohol carbon, was previously resolved arbitrarily by
    BFS/adjacency insertion order rather than by this rule, e.g.
    'CC(O)(Cc1cccc(F)c1)C(F)(F)Br' wrongly named
    '2-(bromodifluoromethyl)-1-(3-fluorophenyl)propan-2-ol' instead of
    PubChem's '1-bromo-1,1-difluoro-3-(3-fluorophenyl)-2-methylpropan-2-ol',
    which absorbs the halogen-bearing carbon into the chain instead,
    citing 5 prefix substituents instead of 2). Every candidate two-arm
    combination (across every pair of distinct starting neighbors of
    `required`, and every farthest node tied within each) is enumerated
    and scored the same way, since `required` rarely has more than two or
    three non-excluded neighbors -- this stays cheap.

    `ring_boundary`/`excluded`/`halogens`: same meaning as
    `longest_branched_chain`.

    Returns (chain, branches): `chain` is the winning path as an
    atom-index list, in an arbitrary direction -- the caller tries both
    ways (same as it already does for a plain unbranched chain) to give
    `required` its own lowest locant; `branches` is {position ->
    [branch_root_atom, ...]}, 1-based against this `chain`'s order, same
    shape as `longest_branched_chain`."""
    blocked = set(ring_boundary) | set(excluded) | set(halogens)
    neighbors = [n for n in graph[required] if n not in blocked]

    groups = []
    for start in neighbors:
        # `required` itself must stay off-limits here -- otherwise this
        # arm's search loops back through it into the *other* arm(s),
        # corrupting the two-arm split with a self-crossing path.
        arm_blocked = blocked | {required}
        dist = {start: 0}
        parent = {start: None}
        queue = [start]
        while queue:
            next_queue = []
            for node in queue:
                for neighbor in graph[node]:
                    if neighbor in arm_blocked or neighbor in dist:
                        continue
                    dist[neighbor] = dist[node] + 1
                    parent[neighbor] = node
                    next_queue.append(neighbor)
            queue = next_queue
        far = max(dist.values())
        groups.append([
            [required] + path_between(parent, start, node)
            for node, d in dist.items()
            if d == far
        ])

    def candidate_chains():
        if not groups:
            yield [required]
        elif len(groups) == 1:
            for arm in groups[0]:
                yield list(reversed(arm))
        else:
            for i in range(len(groups)):
                for j in range(len(groups)):
                    if i == j:
                        continue
                    for arm1 in groups[i]:
                        for arm2 in groups[j]:
                            yield list(reversed(arm1)) + arm2[1:]

    def branches_for(chain):
        chain_set = set(chain)
        branches = {}
        for position, atom in enumerate(chain, start=1):
            # Any neighbor not on the two chosen arms -- including a
            # third+ arm off `required` itself, e.g. a ketone carbon with
            # more than two carbon substituents (shouldn't normally
            # arise, but falls through safely here rather than being
            # silently dropped) -- is an ordinary branch, named via
            # `name_branch` like any other.
            roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
            if roots:
                branches[position] = roots
        return branches

    best_chain, best_branches, best_key = None, None, None
    for chain in candidate_chains():
        branches = branches_for(chain)
        substituent_count = sum(len(roots) for roots in branches.values())
        locants = lowest_locant_set(pos for pos, roots in branches.items() for _ in roots)
        key = (-len(chain), -substituent_count, locants)
        if best_key is None or key < best_key:
            best_key, best_chain, best_branches = key, chain, branches
    return best_chain, best_branches


def halogen_substituents(mol):
    """{atom_idx -> substituent prefix name} for every halogen atom in `mol`
    (P-35.2.1). Passed down into `name_branch` so it can name a halogen leaf
    and exclude halogens from a compound substituent's own internal chain
    search, the same way `carbon_adjacency` does for a parent hydride."""
    return {
        atom.GetIdx(): HALOGEN_PREFIXES[atom.GetAtomicNum()]
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() in HALOGEN_PREFIXES
    }


def bfs(graph, start):
    dist = {start: 0}
    parent = {start: None}
    queue = [start]
    while queue:
        next_queue = []
        for node in queue:
            for neighbor in graph[node]:
                if neighbor not in dist:
                    dist[neighbor] = dist[node] + 1
                    parent[neighbor] = node
                    next_queue.append(neighbor)
        queue = next_queue
    return dist, parent


def component_subgraph(graph, start):
    """`graph` restricted to the connected component containing `start`
    (e.g. one acyl branch of a symmetric multi-fragment structure once the
    bridging atom(s) are excluded from `graph` itself)."""
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def path_between(parent, start, end):
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    return list(reversed(path))


def ring_cycle(graph, ring_atoms):
    """Order a monocyclic ring's atoms into a single walk around the ring,
    starting from `ring_atoms[0]` (an arbitrary RDKit ring-atom listing
    order) -- both traversal directions are still tried by the caller when
    picking the winning numbering."""
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def linear_branch(graph, root, coming_from):
    """Walk a branch outward; return its atom count, or None if it forks
    (a "compound" substituent, P-29.4, not yet supported)."""
    length = 1
    previous, current = coming_from, root
    while True:
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) == 0:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]
        length += 1


def unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None
    if the chain branches, rings, or leaves carbon at any point."""
    length = 0
    previous = exclude_idx
    current = root_idx
    while True:
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return None
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous]
        length += 1
        if not neighbors:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]


def lowest_locant_set(locants):
    return tuple(sorted(locants))


def multiplied_word(count, base):
    """Multiplying-prefix word for `count` occurrences of a suffix like 'ol'/
    'one'/'al' (P-14.2.1): omitted for zero, bare for one, else prefixed with
    the basic numerical term ('di', 'tri', ...).

    P-16.3.3: a multiplying prefix's terminal 'a' (tetra, penta, hexa, ...)
    is elided before a suffix beginning with 'a' or 'o' -- e.g. 'tetra' +
    'ol' -> 'tetrol', not 'tetraol' (PubChem CID 8998, confirming
    "butane-1,2,3,4-tetrol"); 'tetra' + 'amine' -> 'tetramine' (CID
    6395580). 'di'/'tri' never end in 'a', so they're never affected."""
    if count == 0:
        return ""
    if count == 1:
        return base
    prefix = numerical_term(count)
    if prefix.endswith("a") and base[:1] in "ao":
        prefix = prefix[:-1]
    return prefix + base


def elides_before(word: str) -> bool:
    """True if a preceding word ending in 'e' (e.g. an 'ene' locant-suffix
    segment assembled by `multiplied_word`) should have that 'e' elided
    right before `word`. P-31.1.1.1's 'dien-...-yne' elision is specific to
    a following yne/diyne/triyne/... segment -- it always applies there
    regardless of the multiplying prefix, since 'yne' itself is 'y'-initial
    (found via real-data testing: 'deca-1,2,3-trien-5,7,9-triyne' and
    'methyl 18-bromooctadeca-9,17-dien-5,7-diynoate' were wrongly left
    unelided by the ~30 suffix modules that had inlined this check against
    the final word's own first letter). For every other following word
    (oic, diazonium, oate, ...) elision follows that word's own literal
    first letter, not its underlying un-prefixed form -- 'but-2-enedioic
    acid' and 'pent-4-ene-1-diazonium' both stay unelided even though
    'dioic'/'diazonium' start with a consonant only because of their own
    multiplying prefix."""
    if word.endswith("yne"):
        return True
    return word[:1] in "aeiouy"


def group_substituents(substituents):
    """{position -> [(name, is_compound), ...]} -> {name -> {"locants": [...],
    "compound": bool}}, merging same-named substituents at different
    positions so they can be cited once with a multiplying prefix (P-16)."""
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def longest_chains(graph):
    """Every longest path (by atom count) through an undirected acyclic
    `graph` (P-44.3, the principal chain must be a longest chain candidate)."""
    nodes = list(graph)
    distances = {}
    parents = {}
    for node in nodes:
        dist, parent = bfs(graph, node)
        distances[node] = dist
        parents[node] = parent

    diameter = max(d for dist in distances.values() for d in dist.values())
    chains = []
    seen = set()
    for u in nodes:
        for v, d in distances[u].items():
            if d == diameter and (v, u) not in seen:
                seen.add((u, v))
                chains.append(path_between(parents[u], u, v))
    return chains


def bond_locant(chain, bond_atoms):
    """1-based position along `chain` of the bond between `bond_atoms` (a
    2-tuple of atom indices), or None if that bond doesn't lie on `chain`."""
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def bond_locants(chain, bonds):
    """(ene_locants, yne_locants) for every (a, b, order) bond in `bonds` that
    lies on `chain`, or None if any bond doesn't (a shorter/wrong chain
    candidate)."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = bond_locant(chain, (a, b))
        if locant is None:
            return None
        (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return ene, yne


def ring_bond_locant(position_of, bond_atoms, ring_size):
    """1-based ring locant of the lower-numbered atom of a ring bond, under
    `position_of` (atom index -> 1-based position), wrapping so the bond
    between positions 1 and `ring_size` reports as `ring_size`."""
    pa, pb = position_of[bond_atoms[0]], position_of[bond_atoms[1]]
    return ring_size if {pa, pb} == {1, ring_size} else min(pa, pb)


def ring_bond_locants(position_of, bonds, ring_size):
    """(ene_locants, yne_locants), both sorted, for every ring C=C/C#C bond
    under this ring numbering."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = ring_bond_locant(position_of, (a, b), ring_size)
        (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return sorted(ene), sorted(yne)


def specified_stereocenters(mol):
    """Scan `mol` for stereo elements (RDKit's `Chem.FindPotentialStereo`).
    If there are none at all, or every one present is left unspecified
    (no `@`/`@@`/E-Z bond marker anywhere in the input), return None --
    the caller should proceed exactly as if stereochemistry weren't a
    factor. This matches both this project's existing, long-standing
    behavior and the real IUPAC/PubChem convention for a name that doesn't
    specify configuration at all (e.g. this project's own pre-existing,
    PubChem-verified '2-fluorobutan-1-ol', whose C2 is a genuine but
    undrawn stereocenter) -- it is deliberately NOT treated as a new
    rejection case, unlike a *partially* specified molecule (see below).

    If one or more stereo elements are *specified*, no unspecified one
    alongside them, and every specified one is a tetrahedral atom
    stereocenter (no double-bond E/Z element), return a list of
    (atom_idx, "R" or "S") pairs, one per specified stereocenter, via
    `rdCIPLabeler` -- in the same order `Chem.FindPotentialStereo` reports
    them, not yet locant-sorted (the caller only learns each atom's locant
    once its own chain/ring numbering is fixed; see P-91.3's ascending-
    locant citation order in the module docstring). Otherwise -- a
    specified stereocenter mixed with an unspecified one, or any
    double-bond E/Z stereo element -- raise `UnsupportedStructure`
    explicitly (P-93 double-bond stereo combined with P-92 tetrahedral
    stereo, and any partially-specified molecule, are both out of
    scope)."""
    elements = Chem.FindPotentialStereo(mol)
    specified = [e for e in elements if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    if len(specified) != len(elements) or any(e.type != Chem.StereoType.Atom_Tetrahedral for e in specified):
        raise UnsupportedStructure(
            "stereochemistry beyond one or more specified tetrahedral "
            "stereocenters (with no unspecified one alongside them) is not "
            "supported yet (a specified stereocenter mixed with an "
            "unspecified one, or any C=C/C#N double-bond E/Z stereo -- see "
            "P-92/P-93)"
        )
    rdCIPLabeler.AssignCIPLabels(mol)
    labels = []
    for element in specified:
        atom_idx = element.centeredOn
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.HasProp("_CIPCode"):
            raise UnsupportedStructure(
                "could not determine a CIP R/S label for this stereocenter"
            )
        code = atom.GetProp("_CIPCode")
        if code not in ("R", "S"):
            raise UnsupportedStructure(
                "a pseudoasymmetric stereocenter (lowercase 'r'/'s') is not "
                "supported yet -- only uppercase R/S stereocenters are in "
                "scope"
            )
        labels.append((atom_idx, code))
    return labels


def specified_double_bond_stereo(mol):
    """`Bond_Double` analogue of `specified_stereocenters` above (see its
    docstring for the same unspecified-vs-rejected reasoning). None if
    there's no *specified* double-bond stereo element at all -- covering
    both a non-stereogenic double bond (e.g. `C=C(C)C`, where
    `FindPotentialStereo` doesn't report an element at all) and one left
    unspecified in the input (plain `C=C`, no `/`/`\\`) -- so the caller
    proceeds exactly as before (no E/Z prefix) in either case.

    If one or more stereo elements are *specified*, no unspecified one
    alongside them, and every specified one is a `Bond_Double` element (no
    tetrahedral atom stereocenter), return a list of
    (bond_idx, "E" or "Z") pairs via `rdCIPLabeler` -- in the same order
    `Chem.FindPotentialStereo` reports them, not yet locant-sorted (the
    caller only learns each bond's locant once its own chain numbering is
    fixed). Otherwise -- a specified element mixed with an unspecified one,
    or a tetrahedral stereocenter -- raise `UnsupportedStructure`
    explicitly."""
    elements = Chem.FindPotentialStereo(mol)
    specified = [e for e in elements if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    if len(specified) != len(elements) or any(e.type != Chem.StereoType.Bond_Double for e in specified):
        raise UnsupportedStructure(
            "stereochemistry beyond one or more specified C=C double-bond "
            "E/Z elements (with no unspecified one alongside them) is not "
            "supported yet (a specified element mixed with an unspecified "
            "one, or a tetrahedral stereocenter -- see P-92/P-93)"
        )
    rdCIPLabeler.AssignCIPLabels(mol)
    labels = []
    for element in specified:
        bond_idx = element.centeredOn
        bond = mol.GetBondWithIdx(bond_idx)
        if not bond.HasProp("_CIPCode"):
            raise UnsupportedStructure(
                "could not determine a CIP E/Z label for this double bond"
            )
        labels.append((bond_idx, bond.GetProp("_CIPCode")))
    return labels


def specified_stereo_elements(mol):
    """Like `specified_stereocenters`/`specified_double_bond_stereo`, but
    allows a specified tetrahedral stereocenter and a specified C=C
    double-bond E/Z element to coexist in the same molecule: P-91.3's own
    worked example, '(2Z,5R,7E)-nona-2,7-dien-5-ol (PIN)', cites both kinds
    together in one locant-ascending group, confirmed directly from the
    primary source text.

    Returns None if there is no specified stereo element at all (same
    unspecified-vs-rejected policy as the two functions above). Otherwise
    returns a list of ("atom" or "bond", atom_or_bond_idx, "R"/"S"/"E"/"Z")
    triples, one per specified element, not yet locant-sorted -- the
    caller only learns each element's locant once its own chain numbering
    is fixed. Raises `UnsupportedStructure` for an unspecified element
    mixed in, a pseudoasymmetric (lowercase r/s) stereocenter, or a CIP
    label RDKit couldn't determine."""
    elements = Chem.FindPotentialStereo(mol)
    specified = [e for e in elements if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    if len(specified) != len(elements):
        raise UnsupportedStructure(
            "a specified stereo element alongside an unspecified one is "
            "not supported yet (see P-92/P-93)"
        )
    rdCIPLabeler.AssignCIPLabels(mol)
    labels = []
    for element in specified:
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            atom = mol.GetAtomWithIdx(element.centeredOn)
            if not atom.HasProp("_CIPCode"):
                raise UnsupportedStructure(
                    "could not determine a CIP R/S label for this stereocenter"
                )
            code = atom.GetProp("_CIPCode")
            if code not in ("R", "S"):
                raise UnsupportedStructure(
                    "a pseudoasymmetric stereocenter (lowercase 'r'/'s') is "
                    "not supported yet -- only uppercase R/S stereocenters "
                    "are in scope"
                )
            labels.append(("atom", element.centeredOn, code))
        elif element.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(element.centeredOn)
            if not bond.HasProp("_CIPCode"):
                raise UnsupportedStructure(
                    "could not determine a CIP E/Z label for this double bond"
                )
            labels.append(("bond", element.centeredOn, bond.GetProp("_CIPCode")))
        else:
            raise UnsupportedStructure(
                "stereochemistry beyond a tetrahedral R/S stereocenter or a "
                "C=C double-bond E/Z element is not supported yet (see "
                "P-92/P-93)"
            )
    return labels
