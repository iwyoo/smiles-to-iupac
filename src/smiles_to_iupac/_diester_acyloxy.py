"""Naming of a polyester on a plain unbranched polyol chain, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-65.6.3.3.3.1: when all the acyl groups are **identical** (two or
  more), the PIN cites the polyol backbone as a multivalent organyl
  group ('...diyl'/'...triyl'/'...tetrayl'/...) followed by the
  multiplied anion name -- the chapter's own worked examples are
  'ethane-1,2-diyl diacetate (PIN)' and 'propane-1,2,3-triyl triacetate
  (PIN)' (`tmp/bluebook/P6.txt` ~line 7118/7256). This project's
  systematic-name convention (see `_ester.py`'s own 'methyl methanoate'/
  'methyl ethanoate', not 'methyl formate'/'methyl acetate') carries over
  here too, so this module produces 'ethane-1,2-diyl diethanoate'/
  'propane-1,2,3-triyl triethanoate'-style names instead of PubChem's
  retained-name 'diacetate'/'triacetate'.
- P-65.6.3.3.3.2 (method 1): when the acyl groups **differ** (two or
  more), the PIN instead cites each distinct acid as its own anion name,
  in alphanumerical order, each with its own locant set and multiplying
  prefix (`di`/`tri`/...) if it covers more than one position -- e.g.
  'propane-1,2,3-triyl 1,2-diacetate 3-propanoate (PIN)'. For exactly two
  esters, no acid locant is cited at all: the backbone's own two-fold
  symmetry (either numbering direction describes the identical molecule)
  makes one redundant, confirmed directly from the Blue Book's own
  'methylene acetate formate (PIN)' and '1,4-phenylene acetate
  dichloroacetate (PIN)' examples, neither of which cites one. Three or
  more positions lose that symmetry (an interior position is
  constitutionally distinct from the ends), so every group's locant(s)
  are always cited there. Method 2 (citing one ester as the suffix parent
  and the rest as '<acid>yloxy' prefixes) is sanctioned general
  nomenclature only, not implemented here -- this module always returns
  the method-1 PIN for the differing-acid case.
- Scope, deliberately narrow: two or more ester groups, all acyl chains
  plain (unbranched, saturated, no halogens/stereocenters), and every
  ester's alcohol-side carbon lying on a single plain carbon chain --
  either fully unbranched, or, for the exactly-two-identical-esters case
  only, carrying up to one plain saturated alkyl substituent (no ester,
  no ring, no unsaturation, no heteroatom) hanging off a non-terminal
  backbone carbon (e.g. 'propane-1,2-diyl diethanoate',
  '2-methylpropane-1,3-diyl diethanoate'). Three or more identical
  esters with a branch, and a differing-acid diester with a branch, are
  both still out of scope (raise `UnsupportedStructure`) -- the branch's
  own interaction with per-position acid citation (P-65.6.3.3.3.2 method
  1) isn't verified here. A backbone with a branch
  point/quaternary carbon (e.g. pentaerythritol's tetrahedral center) is
  its own separate, not-yet-scoped shape -- not attempted here even for
  identical acids, since it isn't a single chain at all. A branched
  backbone, a branched or unsaturated/halogenated acyl chain, and any
  ring are all out of scope and raise `UnsupportedStructure`.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, carbon_adjacency, non_single_bonds, ordered_chain, path_between
from ._numerals import alkane_name, alkyl_name, multiplying_prefix

_ALLOWED_ATOMIC_NUMS = {6, 8}


def _find_ester_carbons(mol):
    """Every carbon shaped like an ester acyl carbon (a carbonyl oxygen plus
    a second, carbon-bonded ester oxygen), as (acyl_carbon, carbonyl_oxygen,
    ester_oxygen, alcohol_carbon) tuples."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_oxygens = [
            o
            for o in oxygens
            if o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        ]
        if len(carbonyls) == 1 and len(ester_oxygens) == 1:
            alcohol_carbon = next(n for n in ester_oxygens[0].GetNeighbors() if n.GetIdx() != atom.GetIdx())
            matches.append((atom, carbonyls[0], ester_oxygens[0], alcohol_carbon))
    return matches


def has_diester_shape(mol) -> bool:
    return len(_find_ester_carbons(mol)) >= 2


def _tree_path(sub, start, end):
    parent = {start: None}
    queue = [start]
    while queue:
        node = queue.pop(0)
        for neighbor in sub[node]:
            if neighbor not in parent:
                parent[neighbor] = node
                queue.append(neighbor)
    return path_between(parent, start, end)


def _find_full_chain(mol, graph, matches):
    """(chain, branch, attachment_atom) for the polyol backbone carrying
    every one of `matches`' alcohol-side carbons: `chain` is the main
    carbon path (as an atom-index list), `branch` is either `[]` (a fully
    unbranched backbone) or the plain alkyl substituent's own atom-index
    list (its own root first, not including `attachment_atom`), and
    `attachment_atom` (`None` when `branch` is empty) is the specific
    `chain` atom the branch hangs off -- or None (just the one value, not
    a triple) if no such shape exists at all.

    Found via `carbon_adjacency` (ester oxygens are never part of this
    graph at all, so they can never masquerade as a branch) rather than a
    single-direction walk from one ester: the backbone's own connected
    carbon-only component is a tree with either no degree-3+ node (a
    plain chain -- its own two leaves are exactly the chain's two ends,
    including any plain terminal carbon beyond the outermost ester, e.g.
    propanediol diacetate's own terminal methyl) or exactly one (the
    branch point), in which case that tree always has exactly three
    leaves, and exactly one of the three ways to pick two of them as the
    main chain's own ends leaves every ester on the chain and none on the
    leftover leaf's own path (the branch) -- more than one degree-3+ node,
    or a degree-4+ node (two or more branches off one atom), is out of
    scope and returns None here."""
    alcohol_atoms = {m[3].GetIdx() for m in matches}
    carbon_graph = carbon_adjacency(mol)

    component = set()
    stack = [next(iter(alcohol_atoms))]
    while stack:
        node = stack.pop()
        if node in component:
            continue
        component.add(node)
        stack.extend(carbon_graph[node])
    if not alcohol_atoms <= component:
        return None

    sub = {node: [n for n in carbon_graph[node] if n in component] for node in component}
    degree = {node: len(neighbors) for node, neighbors in sub.items()}
    if any(d >= 4 for d in degree.values()):
        return None
    branch_points = [node for node, d in degree.items() if d >= 3]
    if len(branch_points) > 1:
        return None
    # Sorted by atom index so a fully symmetric backbone (every position
    # substituted, both numbering directions equally valid per P-14.4)
    # gets a single deterministic direction here -- `name_diester_acyloxy`'s
    # own candidate loop still tries the reverse too, for the ordinary case
    # where the two directions genuinely differ in their locant set.
    leaves = sorted(node for node, d in degree.items() if d <= 1)

    if not branch_points:
        if len(component) == 1:
            return list(component), [], None
        if len(leaves) != 2:
            return None
        chain = _tree_path(sub, leaves[0], leaves[1])
        if set(chain) != component or not alcohol_atoms <= set(chain):
            return None
        return chain, [], None

    (branch_point,) = branch_points
    if len(leaves) != 3:
        return None
    for branch_leaf in leaves:
        chain_leaves = [leaf for leaf in leaves if leaf != branch_leaf]
        chain = _tree_path(sub, chain_leaves[0], chain_leaves[1])
        branch_path = _tree_path(sub, branch_point, branch_leaf)[1:]
        if not alcohol_atoms <= set(chain):
            continue
        if alcohol_atoms & set(branch_path):
            continue
        if set(chain) | set(branch_path) != component:
            continue
        return chain, branch_path, branch_point
    return None


def _acyl_chain_length(mol, acyl_carbon, ester_oxygen_idx, carbonyl_oxygen_idx):
    graph = adjacency(mol)
    chain = ordered_chain(graph, acyl_carbon.GetIdx(), ester_oxygen_idx, {carbonyl_oxygen_idx})
    if chain is None or any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in chain):
        raise UnsupportedStructure(
            "a branched, unsaturated, or halogen-bearing acyl chain in a "
            "diester is not supported yet (see module docstring)"
        )
    return len(chain)


def _acid_stem(chain_length: int) -> str:
    return alkane_name(chain_length)[:-1]


def name_diester_acyloxy(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the esters' own oxygens are not "
                "supported yet for a polyester (see module docstring)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic ring is out of scope for this module")
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")

    matches = _find_ester_carbons(mol)
    if len(matches) < 2:
        raise UnsupportedStructure("this module only handles two or more ester groups")

    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    if total_oxygens != 2 * len(matches):
        raise UnsupportedStructure(
            "an oxygen outside the esters' own carbonyl/ester pairs is "
            "out of scope for this module"
        )
    carbonyl_bonds = {frozenset((acyl.GetIdx(), carbonyl.GetIdx())) for acyl, carbonyl, _, _ in matches}
    if any(frozenset((a, b)) not in carbonyl_bonds for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure("unsaturation alongside a polyester's ester groups is not supported yet")

    graph = adjacency(mol)
    found = _find_full_chain(mol, graph, matches)
    if found is None:
        raise UnsupportedStructure(
            "the esters must share a single plain carbon chain, with at "
            "most one plain alkyl branch, and no other substituent (see "
            "module docstring)"
        )
    backbone, branch, attachment_atom = found

    lengths = [_acyl_chain_length(mol, acyl, ester_o.GetIdx(), carbonyl.GetIdx()) for acyl, carbonyl, ester_o, _ in matches]

    if branch and (len(matches) != 2 or len(set(lengths)) != 1):
        raise UnsupportedStructure(
            "a plain alkyl branch on the backbone is only supported for "
            "exactly two identical-acid esters (see module docstring)"
        )

    alcohol_atoms = {alcohol.GetIdx() for _, _, _, alcohol in matches}
    best_key = None
    best_position_of = None
    for candidate in (backbone, list(reversed(backbone))):
        position_of = {idx: i + 1 for i, idx in enumerate(candidate)}
        locants = sorted(position_of[idx] for idx in alcohol_atoms)
        key = (locants, position_of[attachment_atom]) if branch else (locants,)
        if best_key is None or key < best_key:
            best_key, best_position_of = key, position_of
    best_locants = best_key[0]

    yl_prefix = multiplying_prefix(len(matches))
    locant_str = ",".join(str(loc) for loc in best_locants)
    group_name = f"{alkane_name(len(backbone))}-{locant_str}-{yl_prefix}yl"
    if branch:
        branch_locant = best_position_of[attachment_atom]
        group_name = f"{branch_locant}-{alkyl_name(len(branch))}{group_name}"

    if len(set(lengths)) == 1:
        acid_name = f"{yl_prefix}{_acid_stem(lengths[0])}oate"
        return f"{group_name} {acid_name}"

    # P-65.6.3.3.3.2 method 1: differing acyl groups are cited as separate
    # anion names, each with its own locant set and multiplying prefix (if
    # repeated), in alphanumerical order of the acid name -- the PIN, in
    # place of the acyloxy-prefix method 2 this module used to return here
    # (general nomenclature only, still `_diester_acyloxy.py`'s territory
    # for reference in the module docstring, not this function's output).
    groups = {}
    for length, (_, _, _, alcohol) in zip(lengths, matches):
        groups.setdefault(length, []).append(best_position_of[alcohol.GetIdx()])

    # For exactly two esters, the backbone's own two-fold symmetry (either
    # numbering direction describes the identical molecule) means no
    # locant is needed to know which acid sits where -- confirmed against
    # the Blue Book's own 'methylene acetate formate (PIN)' and
    # '1,4-phenylene acetate dichloroacetate (PIN)' examples, neither of
    # which cites an acid locant. Three or more positions lose that
    # symmetry (an interior position is constitutionally distinct from
    # the ends), so every group's locant(s) must be cited there, per the
    # 'propane-1,2,3-triyl 1,2-diacetate 3-propanoate (PIN)' example.
    cite_locants = len(matches) != 2

    acid_parts = []
    for length in sorted(groups, key=lambda length: _acid_stem(length)):
        locants = sorted(groups[length])
        count = len(locants)
        prefix = multiplying_prefix(count) if count > 1 else ""
        acid_name = f"{prefix}{_acid_stem(length)}oate"
        if cite_locants:
            loc_str = ",".join(str(loc) for loc in locants)
            acid_name = f"{loc_str}-{acid_name}"
        acid_parts.append(acid_name)

    return f"{group_name} " + " ".join(acid_parts)
