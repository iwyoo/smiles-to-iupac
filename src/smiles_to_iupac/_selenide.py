"""Naming of selenides (the 'selanyl' substituent prefix, R-Se-R'),
restricted to two acyclic saturated hydrocarbon chains hung off a single
selenide selenium, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-63.2.1: a selenide, like a sulfide (`_sulfide.py`) or ether
  (`_ether.py`), has no principal-characteristic-group suffix, so
  R-Se-R' is named substitutively with one side (RH) as parent hydride and
  the other (R'-Se-) as a substituent prefix on it, with parent choice
  following the same P-44.3 skeletal-atom-count rule. This module mirrors
  `_sulfide.py`'s structure directly, selenium in place of sulfur and the
  'selanyl' prefix in place of 'sulfanyl' -- confirmed via PubChem PUG
  REST: CID 11648 (`C[Se]C`) -> "methylselanylmethane", CID 61173
  (`CC[Se]CC`) -> "ethylselanylethane", CID 12248622 (`C[Se]CC`) ->
  "methylselanylethane" (the shorter methyl side becomes the substituent,
  the longer ethyl side the parent), CID 15932888 (`CCC[Se]C`) ->
  "1-methylselanylpropane" (a 3-carbon parent needs the locant -- the same
  citation rule `_sulfide.py` already gets for free from the shared
  `name_from_carbon_graph` helper, confirmed directly against that
  module's own identical `CCCSC` -> "1-methylsulfanylpropane").

Scope and out-of-scope structures are identical to `_sulfide.py`, selenium
in place of sulfur -- see that module's docstring, including the
P-63.2.2.1.1 enclosure pattern for a branched R' substituent. Still out of
scope: any unsaturation or ring, and any heteroatom other than the single
selenide selenium (in particular a diselenide Se-Se, or an oxidized
selenium -- selenoxide/selenone -- are separate functional groups, not in
scope here).
"""

from ._acyclic import longest_chain_length, name_from_carbon_graph, winning_chain_with_key
from ._common import UnsupportedStructure, adjacency, bfs, carbon_adjacency, non_single_bonds
from ._substituents import name_branch

_SELENIUM = 34


def _selanyl_prefix(name: str) -> str:
    return name + "selanyl"


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_selenide_shape(mol) -> bool:
    """True iff `mol` has exactly one selenium atom, singly bonded to two
    carbons (a plain selenide -Se-, P-63.2.1) -- the shape this module
    accepts."""
    seleniums = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _SELENIUM]
    if len(seleniums) != 1:
        return False
    (selenium,) = seleniums
    return selenium.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in selenium.GetNeighbors())


def name_selenide(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, _SELENIUM):
            raise UnsupportedStructure(
                "heteroatoms other than the selenide selenium are not "
                "supported yet (P-63.2.1 is restricted to a plain -Se- "
                "selenide here)"
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
            "alkenes/alkynes; not yet combined with a selenide here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    (selenium,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _SELENIUM)
    selenium_idx = selenium.GetIdx()
    n1, n2 = (n.GetIdx() for n in selenium.GetNeighbors())

    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    size1 = len(_component_subgraph(carbon_graph, n1))
    size2 = len(_component_subgraph(carbon_graph, n2))

    if size1 == size2:
        graph_a = _component_subgraph(carbon_graph, n1)
        graph_b = _component_subgraph(carbon_graph, n2)
        len_a, len_b = longest_chain_length(graph_a), longest_chain_length(graph_b)
        if len_a > len_b:
            parent_root, sub_root = n1, n2
        elif len_b > len_a:
            parent_root, sub_root = n2, n1
        else:
            name_a, compound_a = name_branch(full_graph, n1, selenium_idx, {})
            name_b, compound_b = name_branch(full_graph, n2, selenium_idx, {})
            sub_from_a = f"({name_a})" if compound_a else name_a
            sub_from_b = f"({name_b})" if compound_b else name_b
            key_a, _, _ = winning_chain_with_key(full_graph, graph_a, {selenium_idx: _selanyl_prefix(sub_from_b)})
            key_b, _, _ = winning_chain_with_key(full_graph, graph_b, {selenium_idx: _selanyl_prefix(sub_from_a)})
            parent_root, sub_root = (n1, n2) if key_a <= key_b else (n2, n1)
    elif size1 > size2:
        parent_root, sub_root = n1, n2
    else:
        parent_root, sub_root = n2, n1

    sub_name, sub_compound = name_branch(full_graph, sub_root, selenium_idx, {})
    if sub_compound:
        sub_name = f"({sub_name})"

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {selenium_idx: _selanyl_prefix(sub_name)}
    return name_from_carbon_graph(full_graph, parent_carbon_graph, terminals)
