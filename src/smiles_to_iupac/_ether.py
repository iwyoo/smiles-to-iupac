"""Naming of ethers (the 'oxy' substituent prefix, R-O-R'), restricted to
two acyclic saturated hydrocarbon chains hung off a single ether oxygen, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-63.2.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  ethers have no principal-characteristic-group suffix, so R-O-R' is named
  substitutively with one side (RH) as parent hydride and the other
  (R'-O-) as a substituent prefix on it. Choice of parent follows the
  general parent-hydride seniority rules (P-44), which for two acyclic
  saturated chains reduces to P-44.3: the side with the greater number of
  skeletal (carbon) atoms is the parent.
- P-63.2.2.1.1: the R'-O- substituent prefix is formed by adding 'oxy' to
  the substituent-group name for R'. 'methoxy', 'ethoxy', 'propoxy', and
  'butoxy' are the retained contracted forms for the four shortest
  unbranched chains; any other unbranched chain uses the full alkyl name
  plus 'oxy' unchanged (e.g. 'pentyloxy').
- A branched (compound) R' encloses only R' in parentheses, with 'oxy'
  outside them (e.g. '(butan-2-yl)oxy', per the Blue Book's own worked
  example) — a different enclosure pattern than every other compound
  substituent in this project, which parenthesizes its whole name
  including any trailing text (see `_substituents.py`). R (the parent
  side) may still be branched, since its own substituents are named by
  the ordinary `_acyclic.py` machinery, unaffected by this restriction.
- P-91.3/P-92 (`tasks/ether-stereocenter-naming.md`): a molecule with one
  or more *specified* tetrahedral stereocenters on the parent (R) chain
  gets a "(<locant><R/S>,...)-" prefix, ascending locant order, e.g.
  '(2S)-2-ethoxybutane' -- same mechanism as `_acetal.py`, using
  `_acyclic.py`'s `winning_chain_from_carbon_graph` for the parent
  chain's locant lookup. A stereocenter on the alkoxy (R') substituent
  branch remains out of scope (raises `UnsupportedStructure`).

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one oxygen, or an oxygen not shaped like a plain ether (degree
  != 2, a non-carbon neighbor) — routed to a different module by `core.py`
  before this one is even tried.
- Any unsaturation, any ring, or any heteroatom other than the single ether
  oxygen.
- A tie in skeletal-atom count between the two chains where neither side
  can serve as the (necessarily unbranched) R'-substituent side, i.e. both
  chains are branched.
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

_CONTRACTED_OXY = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy"}


def _oxy_prefix(name: str) -> str:
    return _CONTRACTED_OXY.get(name, name + "oxy")


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_ether_shape(mol) -> bool:
    """True iff `mol` has exactly one oxygen atom, singly bonded to two
    carbons (a plain ether -O-, P-63.2.1) — the shape this module accepts."""
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(oxygens) != 1:
        return False
    (oxygen,) = oxygens
    return oxygen.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in oxygen.GetNeighbors())


def name_ether(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8):
            raise UnsupportedStructure(
                "heteroatoms other than the ether oxygen are not supported "
                "yet (P-63.2.1 is restricted to a plain -O- ether here)"
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
            "alkenes/alkynes; not yet combined with an ether here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    (oxygen,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    oxygen_idx = oxygen.GetIdx()
    n1, n2 = (n.GetIdx() for n in oxygen.GetNeighbors())

    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    size1 = len(_component_subgraph(carbon_graph, n1))
    size2 = len(_component_subgraph(carbon_graph, n2))

    if size1 == size2:
        name_a, compound_a = name_branch(full_graph, n1, oxygen_idx, {})
        name_b, compound_b = name_branch(full_graph, n2, oxygen_idx, {})
        if compound_a and compound_b:
            raise UnsupportedStructure(
                "an ether tied in skeletal-atom count with both sides "
                "branched is not supported yet (see P-63.2.2.1.1's "
                "enclosure interaction, module docstring)"
            )
        parent_root, sub_root = (n2, n1) if compound_a else (n1, n2)
    elif size1 > size2:
        parent_root, sub_root = n1, n2
    else:
        parent_root, sub_root = n2, n1

    sub_name, sub_compound = name_branch(full_graph, sub_root, oxygen_idx, {})
    if sub_compound:
        sub_name = f"({sub_name})"

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {oxygen_idx: _oxy_prefix(sub_name)}
    chain, name = winning_chain_from_carbon_graph(full_graph, parent_carbon_graph, terminals)

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter must lie on the winning parent
    # chain; one on the alkoxy (R') substituent branch is out of scope,
    # mirroring `_acetal.py`'s identical treatment.
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    if any(atom not in position_of for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the "
            "principal chain is not supported yet (see P-92)"
        )
    labels = sorted((position_of[atom], code) for atom, code in stereo)
    prefix = ",".join(f"{locant}{code}" for locant, code in labels)
    return f"({prefix})-{name}"
