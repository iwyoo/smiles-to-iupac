"""Naming of tellurides (the 'tellanyl' substituent prefix, R-Te-R'),
restricted to two acyclic saturated hydrocarbon chains hung off a single
telluride tellurium, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-63.2.1: a telluride, like a selenide (`_selenide.py`) or sulfide
  (`_sulfide.py`), has no principal-characteristic-group suffix, so
  R-Te-R' is named substitutively with one side (RH) as parent hydride
  and the other (R'-Te-) as a substituent prefix on it, with parent
  choice following the same P-44.3 skeletal-atom-count rule. This module
  mirrors `_selenide.py`'s structure directly, tellurium in place of
  selenium and the 'tellanyl' prefix in place of 'selanyl' -- confirmed
  via PubChem PUG REST: CID 68977 (`C[Te]C`) -> "methyltellanylmethane",
  CID 69394 (`CC[Te]CC`) -> "ethyltellanylethane", CID 13981584
  (`C[Te]CC`) -> "methyltellanylethane" (the shorter methyl side becomes
  the substituent, the longer ethyl side the parent), CID 15932889
  (`CCC[Te]C`) -> "1-methyltellanylpropane" (a 3-carbon parent needs the
  locant, same rule as `_selenide.py`). Unlike `_diselenide.py`/
  `_ditelluride.py`'s "di(selanyl|tellanyl)" prefixes, 'tellanyl' doesn't
  begin with a multiplying-prefix-look-alike syllable, so this module can
  use the ordinary `_acyclic.name_from_carbon_graph`/`terminals` shortcut
  unchanged, the same way `_selenide.py`/`_sulfide.py` do.

Scope and out-of-scope structures are identical to `_selenide.py`,
tellurium in place of selenium -- see that module's docstring, including
the P-63.2.2.1.1 enclosure pattern for a branched R' substituent. Still
out of scope: any unsaturation or ring, and any heteroatom other than the
single telluride tellurium (in particular a ditelluride Te-Te, or
oxidized tellurium, are separate functional groups, not in scope here).
"""

from ._acyclic import longest_chain_length, name_from_carbon_graph, winning_chain_with_key
from ._common import UnsupportedStructure, adjacency, bfs, carbon_adjacency, non_single_bonds
from ._substituents import name_branch

_TELLURIUM = 52


def _tellanyl_prefix(name: str) -> str:
    return name + "tellanyl"


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_telluride_shape(mol) -> bool:
    """True iff `mol` has exactly one tellurium atom, singly bonded to
    two carbons (a plain telluride -Te-, P-63.2.1) -- the shape this
    module accepts."""
    telluriums = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _TELLURIUM]
    if len(telluriums) != 1:
        return False
    (tellurium,) = telluriums
    return tellurium.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in tellurium.GetNeighbors())


def name_telluride(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, _TELLURIUM):
            raise UnsupportedStructure(
                "heteroatoms other than the telluride tellurium are not "
                "supported yet (P-63.2.1 is restricted to a plain -Te- "
                "telluride here)"
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
            "alkenes/alkynes; not yet combined with a telluride here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    (tellurium,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _TELLURIUM)
    tellurium_idx = tellurium.GetIdx()
    n1, n2 = (n.GetIdx() for n in tellurium.GetNeighbors())

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
            name_a, compound_a = name_branch(full_graph, n1, tellurium_idx, {})
            name_b, compound_b = name_branch(full_graph, n2, tellurium_idx, {})
            sub_from_a = f"({name_a})" if compound_a else name_a
            sub_from_b = f"({name_b})" if compound_b else name_b
            key_a, _, _ = winning_chain_with_key(full_graph, graph_a, {tellurium_idx: _tellanyl_prefix(sub_from_b)})
            key_b, _, _ = winning_chain_with_key(full_graph, graph_b, {tellurium_idx: _tellanyl_prefix(sub_from_a)})
            parent_root, sub_root = (n1, n2) if key_a <= key_b else (n2, n1)
    elif size1 > size2:
        parent_root, sub_root = n1, n2
    else:
        parent_root, sub_root = n2, n1

    sub_name, sub_compound = name_branch(full_graph, sub_root, tellurium_idx, {})
    if sub_compound:
        sub_name = f"({sub_name})"

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {tellurium_idx: _tellanyl_prefix(sub_name)}
    return name_from_carbon_graph(full_graph, parent_carbon_graph, terminals)
