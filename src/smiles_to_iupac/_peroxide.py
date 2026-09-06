"""Naming of peroxides (the 'peroxy' substituent prefix, R-O-O-R'),
restricted to two acyclic saturated hydrocarbon chains hung off a single
-O-O- bridge, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-63.2.5 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  like an ether (`_ether.py`), a peroxide has no principal-characteristic-
  group suffix, so R-O-O-R' is named substitutively with one side (RH) as
  parent hydride and the other (R'-O-O-) as a substituent prefix on it.
  Choice of parent follows the same P-44.3 rule as `_ether.py`: the side
  with the greater number of skeletal (carbon) atoms is the parent.
  Confirmed via a worked example: dimethyl peroxide's PIN is
  'methylperoxymethane', not the retained functional-class name 'dimethyl
  peroxide'.
- The R'-O-O- substituent prefix is formed by adding 'peroxy' to the
  substituent-group name for R' -- unlike 'oxy' (`_ether.py`'s
  'methoxy'/'ethoxy'/'propoxy'/'butoxy'), 'peroxy' has no retained
  contracted forms for short chains, so it is simply concatenated
  (`methyl` + `peroxy` -> `methylperoxy`).
- A branched (compound) R' has the same enclosure pattern as `_ether.py`'s
  alkoxy prefix (P-63.2.2.1.1's own worked example encloses only R', with
  'peroxy' outside), mirroring `_ether.py` exactly.
- P-91.3/P-92: a molecule with
  one or more *specified* tetrahedral stereocenters on the parent (R)
  chain gets a "(<locant><R/S>,...)-" prefix, ascending locant order --
  same mechanism as `_ether.py`, using `_acyclic.py`'s
  `winning_chain_from_carbon_graph` for the parent chain's locant lookup.
  A stereocenter on the peroxy (R') substituent branch remains out of
  scope (raises `UnsupportedStructure`).
- A tie in total skeletal-atom count between the two chains is resolved
  the same two-step way as `_ether.py` (see that module's docstring for
  the worked examples): first by each side's own longest achievable chain
  length (`longest_chain_length`), and only if that also ties, by
  comparing each side's resulting locant set as parent
  (`winning_chain_with_key`).

Explicitly out of scope (raise `UnsupportedStructure`):
- More than two oxygens, or two oxygens not shaped like a plain -O-O-
  bridge (not bonded to each other, degree != 2, a non-carbon second
  neighbor) -- routed to a different module by `core.py` before this one
  is even tried.
- Any unsaturation, any non-benzene ring, or any heteroatom other than the
  peroxide's own two oxygens.

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_peroxide_chain`), mirroring `_ether.py`'s equivalent:
'peroxy' has no suffix form, so P-44.1.2.2 rule (1) makes the ring the
parent regardless of the other side's chain length. Both a direct
ring-oxygen bond and a chain spacer between the ring and the near
peroxide oxygen are supported, confirmed via PubChem PUG REST for the
underlying structure/connectivity (CID 15817998 'c1ccccc1OOCC' ->
'ethylperoxybenzene', 'c1ccccc1COOCC' -> 'ethylperoxymethylbenzene',
'c1ccccc1OOC(C)C' -> 'propan-2-ylperoxybenzene', 'c1ccccc1COOC(C)C' ->
'propan-2-ylperoxymethylbenzene') -- but, like `_ether.py`, a compound R'
side or a compound chain-spacer branch is still parenthesized here even
where PubChem's own auto-generated name omits the parentheses. Unlike
`_disulfide.py`'s narrower benzene-ring path (direct ring-sulfur bond
only), the chain-spacer case is in scope here: 'peroxy' is an independent
prefix word applied only to R', not a fused whole-fragment prefix like
'disulfanyl', so it doesn't hit the nested-parenthesization mismatch that
excluded disulfide's chain-spacer case.
"""

from ._acyclic import longest_chain_length, winning_chain_from_carbon_graph, winning_chain_with_key
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._substituents import name_branch


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def _peroxide_oxygens(mol):
    """(o1, o2) atoms of a plain -O-O- bridge -- two oxygens, singly bonded
    to each other, each with exactly one other neighbor, a carbon -- or
    None if `mol` doesn't have exactly two oxygens shaped this way."""
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(oxygens) != 2:
        return None
    o1, o2 = oxygens
    if o1.GetDegree() != 2 or o2.GetDegree() != 2:
        return None
    bond = mol.GetBondBetweenAtoms(o1.GetIdx(), o2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return None
    for o in (o1, o2):
        (other,) = (n for n in o.GetNeighbors() if n.GetIdx() not in (o1.GetIdx(), o2.GetIdx()))
        if other.GetAtomicNum() != 6:
            return None
    return o1, o2


def has_peroxide_shape(mol) -> bool:
    return _peroxide_oxygens(mol) is not None


def _validate_peroxide_atoms(mol, ring_atoms=frozenset()):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path (mirrors `_ether.py`'s
    `_validate_ether_atoms`): every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a non-aromatic chain carbon or one of
    the peroxide's own two oxygens."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() not in (6, 8):
            raise UnsupportedStructure(
                "heteroatoms other than the peroxide's own two oxygens are "
                "not supported yet (P-63.2.5 is restricted to a plain "
                "-O-O- peroxide here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and idx not in ring_atoms:
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see P-31 for "
            "alkenes/alkynes; not yet combined with a peroxide here)"
        )


def _name_benzene_ring_peroxide_chain(mol, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since the 'peroxy' prefix has no suffix form
    (module docstring), a single, otherwise-unsubstituted benzene ring is
    always the parent hydride, regardless of the other side's chain
    length. Handles both a direct ring-oxygen bond ('ethylperoxybenzene')
    and a chain spacer between the ring and the near peroxide oxygen
    ('(2-ethylperoxyethyl)benzene'), mirroring
    `_ether._name_benzene_ring_ether_chain` with the single ether oxygen
    generalized to the peroxide's near/far oxygen pair."""
    _validate_peroxide_atoms(mol, ring_atoms)
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter is not supported yet for a "
            "benzene-ring-substituent peroxide"
        )

    graph = adjacency(mol)
    o1, o2 = _peroxide_oxygens(mol)
    o1_idx, o2_idx = o1.GetIdx(), o2.GetIdx()

    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a peroxide is not supported yet"
        )
    ring_atom, chain_root = attachment

    if chain_root in (o1_idx, o2_idx):
        near_o = chain_root
        far_o = o2_idx if near_o == o1_idx else o1_idx
        (r_prime,) = [n for n in graph[far_o] if n != near_o]
        sub_name, sub_compound = name_branch(graph, r_prime, far_o, {}, mol=mol)
        if sub_compound:
            sub_name = f"({sub_name})"
        return f"{sub_name}peroxybenzene"

    blocked_graph = {node: [n for n in neighbors if n not in (o1_idx, o2_idx)] for node, neighbors in graph.items()}
    del blocked_graph[o1_idx]
    del blocked_graph[o2_idx]
    reached, _ = bfs(blocked_graph, ring_atom)
    near_o = o1_idx if any(n in reached for n in graph[o1_idx]) else o2_idx
    far_o = o2_idx if near_o == o1_idx else o1_idx
    (r_prime,) = [n for n in graph[far_o] if n != near_o]

    sub_name, sub_compound = name_branch(graph, r_prime, far_o, {}, mol=mol)
    if sub_compound:
        sub_name = f"({sub_name})"
    peroxy_term = sub_name + "peroxy"
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, {near_o: peroxy_term}, mol=mol)
    if not is_compound:
        return f"{branch_name}benzene"
    if "(" in branch_name:
        return f"[{branch_name}]benzene"
    return f"({branch_name})benzene"


def name_peroxide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_peroxide_chain(mol, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_peroxide_atoms(mol)

    o1, o2 = _peroxide_oxygens(mol)
    o1_idx, o2_idx = o1.GetIdx(), o2.GetIdx()
    n1 = next(n.GetIdx() for n in o1.GetNeighbors() if n.GetIdx() != o2_idx)
    n2 = next(n.GetIdx() for n in o2.GetNeighbors() if n.GetIdx() != o1_idx)

    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    size1 = len(_component_subgraph(carbon_graph, n1))
    size2 = len(_component_subgraph(carbon_graph, n2))

    if size1 == size2:
        graph_a = _component_subgraph(carbon_graph, n1)
        graph_b = _component_subgraph(carbon_graph, n2)
        len_a, len_b = longest_chain_length(graph_a), longest_chain_length(graph_b)
        if len_a > len_b:
            parent_root, parent_oxygen, sub_root, sub_oxygen = n1, o1_idx, n2, o2_idx
        elif len_b > len_a:
            parent_root, parent_oxygen, sub_root, sub_oxygen = n2, o2_idx, n1, o1_idx
        else:
            name_a, compound_a = name_branch(full_graph, n1, o1_idx, {}, mol=mol)
            name_b, compound_b = name_branch(full_graph, n2, o2_idx, {}, mol=mol)
            sub_from_a = f"({name_a})" if compound_a else name_a
            sub_from_b = f"({name_b})" if compound_b else name_b
            key_a, _, _ = winning_chain_with_key(full_graph, graph_a, {o1_idx: sub_from_b + "peroxy"}, mol=mol)
            key_b, _, _ = winning_chain_with_key(full_graph, graph_b, {o2_idx: sub_from_a + "peroxy"}, mol=mol)
            parent_root, parent_oxygen, sub_root, sub_oxygen = (
                (n1, o1_idx, n2, o2_idx) if key_a <= key_b else (n2, o2_idx, n1, o1_idx)
            )
    elif size1 > size2:
        parent_root, parent_oxygen, sub_root, sub_oxygen = n1, o1_idx, n2, o2_idx
    else:
        parent_root, parent_oxygen, sub_root, sub_oxygen = n2, o2_idx, n1, o1_idx

    sub_name, sub_compound = name_branch(full_graph, sub_root, sub_oxygen, {}, mol=mol)
    if sub_compound:
        sub_name = f"({sub_name})"

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {parent_oxygen: sub_name + "peroxy"}
    chain, name = winning_chain_from_carbon_graph(full_graph, parent_carbon_graph, terminals, mol=mol)

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter must lie on the winning parent
    # chain; one on the peroxy (R') substituent branch is out of scope,
    # mirroring `_ether.py`'s identical treatment.
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    if any(atom not in position_of for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the "
            "principal chain is not supported yet (see P-92)"
        )
    labels = sorted((position_of[atom], code) for atom, code in stereo)
    prefix = ",".join(f"{locant}{code}" for locant, code in labels)
    return f"({prefix})-{name}"
