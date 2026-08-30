"""Naming of sulfides (the 'sulfanyl' substituent prefix, R-S-R'), restricted
to two acyclic saturated hydrocarbon chains hung off a single sulfide
sulfur, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-63.2.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  sulfides (the sulfur analogue of ethers) have no principal-characteristic-
  group suffix, so R-S-R' is named substitutively with one side (RH) as
  parent hydride and the other (R'-S-) as a substituent prefix on it, with
  parent choice following the same P-44.3 skeletal-atom-count rule as
  `_ether.py`.
- P-63.2.2.1.2: the R'-S- substituent prefix is formed by adding 'sulfanyl'
  to the substituent-group name for R'. Unlike 'oxy', 'sulfanyl' has no
  contracted forms for the short chains -- 'methylsulfanyl', not
  'methsulfanyl' or the older, non-PIN 'methylthio'.

Scope and out-of-scope structures are identical to `_ether.py`, sulfur in
place of oxygen -- see that module's docstring; this one mirrors its
structure directly, including the P-63.2.2.1.1 enclosure pattern for a
branched R' substituent. Still out of scope: any unsaturation or ring, and
any heteroatom other than the single sulfide sulfur (in particular a
disulfide S-S, or an oxidized sulfur -- sulfoxide/sulfone -- are separate
functional groups, not in scope here).

- P-91.3/P-92 (`tasks/sulfide-stereocenter-naming.md`): a molecule with
  one or more *specified* tetrahedral stereocenters on the parent (R)
  chain gets a "(<locant><R/S>,...)-" prefix, ascending locant order --
  same mechanism as `_ether.py`, using `_acyclic.py`'s
  `winning_chain_from_carbon_graph` for the parent chain's locant lookup.
  A stereocenter on the sulfanyl (R') substituent branch remains out of
  scope (raises `UnsupportedStructure`).
"""

from ._acyclic import winning_chain_from_carbon_graph
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    non_single_bonds,
    specified_stereocenters,
)
from ._substituents import name_branch


def _sulfanyl_prefix(name: str) -> str:
    return name + "sulfanyl"


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_sulfide_shape(mol) -> bool:
    """True iff `mol` has exactly one sulfur atom, singly bonded to two
    carbons (a plain sulfide -S-, P-63.2.1) -- the shape this module
    accepts."""
    sulfurs = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if len(sulfurs) != 1:
        return False
    (sulfur,) = sulfurs
    return sulfur.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in sulfur.GetNeighbors())


def name_sulfide(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 16):
            raise UnsupportedStructure(
                "heteroatoms other than the sulfide sulfur are not "
                "supported yet (P-63.2.1 is restricted to a plain -S- "
                "sulfide here)"
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
            "alkenes/alkynes; not yet combined with a sulfide here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    (sulfur,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    sulfur_idx = sulfur.GetIdx()
    n1, n2 = (n.GetIdx() for n in sulfur.GetNeighbors())

    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    size1 = len(_component_subgraph(carbon_graph, n1))
    size2 = len(_component_subgraph(carbon_graph, n2))

    if size1 == size2:
        name_a, compound_a = name_branch(full_graph, n1, sulfur_idx, {})
        name_b, compound_b = name_branch(full_graph, n2, sulfur_idx, {})
        if compound_a and compound_b:
            raise UnsupportedStructure(
                "a sulfide tied in skeletal-atom count with both sides "
                "branched is not supported yet (see P-63.2.2.1.1's "
                "enclosure interaction, module docstring)"
            )
        parent_root, sub_root = (n2, n1) if compound_a else (n1, n2)
    elif size1 > size2:
        parent_root, sub_root = n1, n2
    else:
        parent_root, sub_root = n2, n1

    sub_name, sub_compound = name_branch(full_graph, sub_root, sulfur_idx, {})
    if sub_compound:
        sub_name = f"({sub_name})"

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {sulfur_idx: _sulfanyl_prefix(sub_name)}
    chain, name = winning_chain_from_carbon_graph(full_graph, parent_carbon_graph, terminals)

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter must lie on the winning parent
    # chain; one on the sulfanyl (R') substituent branch is out of scope,
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
