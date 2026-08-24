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
- A branched (compound) R' has the same unimplemented enclosure
  interaction as `_ether.py`'s alkoxy prefix (P-63.2.2.1.1's own worked
  example encloses only R', with 'peroxy' outside), so it is out of scope
  here too, mirroring `_ether.py` exactly.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than two oxygens, or two oxygens not shaped like a plain -O-O-
  bridge (not bonded to each other, degree != 2, a non-carbon second
  neighbor) -- routed to a different module by `core.py` before this one
  is even tried.
- Any unsaturation, any ring, or any heteroatom other than the peroxide's
  own two oxygens.
- A tie in skeletal-atom count between the two chains where neither side
  can serve as the (necessarily unbranched) R'-substituent side, i.e. both
  chains are branched.
"""

from ._acyclic import name_from_carbon_graph
from ._common import UnsupportedStructure, adjacency, bfs, carbon_adjacency, non_single_bonds
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


def name_peroxide(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8):
            raise UnsupportedStructure(
                "heteroatoms other than the peroxide's own two oxygens are "
                "not supported yet (P-63.2.5 is restricted to a plain "
                "-O-O- peroxide here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see P-31 for "
            "alkenes/alkynes; not yet combined with a peroxide here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    o1, o2 = _peroxide_oxygens(mol)
    o1_idx, o2_idx = o1.GetIdx(), o2.GetIdx()
    n1 = next(n.GetIdx() for n in o1.GetNeighbors() if n.GetIdx() != o2_idx)
    n2 = next(n.GetIdx() for n in o2.GetNeighbors() if n.GetIdx() != o1_idx)

    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    size1 = len(_component_subgraph(carbon_graph, n1))
    size2 = len(_component_subgraph(carbon_graph, n2))

    if size1 == size2:
        name_a, compound_a = name_branch(full_graph, n1, o1_idx, {})
        name_b, compound_b = name_branch(full_graph, n2, o2_idx, {})
        if compound_a and compound_b:
            raise UnsupportedStructure(
                "a peroxide tied in skeletal-atom count with both sides "
                "branched is not supported yet (see P-63.2.2.1.1's "
                "enclosure interaction, module docstring)"
            )
        parent_root, parent_oxygen, sub_root, sub_oxygen = (
            (n2, o2_idx, n1, o1_idx) if compound_a else (n1, o1_idx, n2, o2_idx)
        )
    elif size1 > size2:
        parent_root, parent_oxygen, sub_root, sub_oxygen = n1, o1_idx, n2, o2_idx
    else:
        parent_root, parent_oxygen, sub_root, sub_oxygen = n2, o2_idx, n1, o1_idx

    sub_name, sub_compound = name_branch(full_graph, sub_root, sub_oxygen, {})
    if sub_compound:
        raise UnsupportedStructure(
            "a branched peroxy substituent's enclosing marks are not "
            "supported yet (see P-63.2.2.1.1, module docstring)"
        )

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {parent_oxygen: sub_name + "peroxy"}
    return name_from_carbon_graph(full_graph, parent_carbon_graph, terminals)
