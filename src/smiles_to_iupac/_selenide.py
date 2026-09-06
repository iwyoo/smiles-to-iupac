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
scope: any unsaturation or non-benzene ring, and any heteroatom other than
the single selenide selenium (in particular a diselenide Se-Se, or an
oxidized selenium -- selenoxide/selenone -- are separate functional
groups, not in scope here).

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_selenide_chain`), mirroring `_sulfide.py`'s
`_name_benzene_ring_sulfide_chain` exactly (selenium in place of sulfur,
'selanyl' in place of 'sulfanyl'): since 'selanyl' has no suffix form,
P-44.1.2.2 rule (1) makes the ring the parent regardless of the other
side's chain length. Both a direct ring-selenium bond and a chain spacer
are supported, confirmed via PubChem PUG REST (CID 140285
'c1ccccc1[Se]CC' -> 'ethylselanylbenzene', CID 12975596
'c1ccccc1C[Se]CC' -> 'ethylselanylmethylbenzene', CID 140894
'c1ccccc1[Se]C(C)C' -> 'propan-2-ylselanylbenzene').
"""

from rdkit import Chem

from ._acyclic import longest_chain_length, name_from_carbon_graph, winning_chain_with_key
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


def _validate_selenide_atoms(mol, ring_atoms=frozenset()):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path, mirrors `_sulfide.py`'s
    `_validate_sulfide_atoms`."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() not in (6, _SELENIUM):
            raise UnsupportedStructure(
                "heteroatoms other than the selenide selenium are not "
                "supported yet (P-63.2.1 is restricted to a plain -Se- "
                "selenide here)"
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
            "alkenes/alkynes; not yet combined with a selenide here)"
        )


def _name_benzene_ring_selenide_chain(mol, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since the 'selanyl' prefix has no suffix form
    (module docstring), a single, otherwise-unsubstituted benzene ring is
    always the parent hydride, regardless of the other side's chain
    length. Handles both a direct ring-selenium bond
    ('ethylselanylbenzene') and a chain spacer between the ring and the
    selenide selenium ('(2-ethylselanylethyl)benzene'), mirroring
    `_sulfide._name_benzene_ring_sulfide_chain` exactly."""
    _validate_selenide_atoms(mol, ring_atoms)
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter is not supported yet for a "
            "benzene-ring-substituent selenide"
        )

    graph = adjacency(mol)
    (selenium_idx,) = (idx for idx in graph if mol.GetAtomWithIdx(idx).GetAtomicNum() == _SELENIUM)

    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a selenide is not supported yet"
        )
    ring_atom, chain_root = attachment

    if chain_root == selenium_idx:
        (r_prime,) = [n for n in graph[selenium_idx] if n != ring_atom]
        sub_name, sub_compound = name_branch(graph, r_prime, selenium_idx, {}, mol=mol)
        if sub_compound:
            sub_name = f"({sub_name})"
        return f"{_selanyl_prefix(sub_name)}benzene"

    blocked_graph = {node: [n for n in neighbors if n != selenium_idx] for node, neighbors in graph.items()}
    del blocked_graph[selenium_idx]
    reached, _ = bfs(blocked_graph, ring_atom)
    (r_prime,) = [n for n in graph[selenium_idx] if n not in reached]

    sub_name, sub_compound = name_branch(graph, r_prime, selenium_idx, {}, mol=mol)
    if sub_compound:
        sub_name = f"({sub_name})"
    selanyl_term = _selanyl_prefix(sub_name)
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, {selenium_idx: selanyl_term}, mol=mol)
    if not is_compound:
        return f"{branch_name}benzene"
    if "(" in branch_name:
        return f"[{branch_name}]benzene"
    return f"({branch_name})benzene"


def name_selenide(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_selenide_chain(mol, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_selenide_atoms(mol)

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
            name_a, compound_a = name_branch(full_graph, n1, selenium_idx, {}, mol=mol)
            name_b, compound_b = name_branch(full_graph, n2, selenium_idx, {}, mol=mol)
            sub_from_a = f"({name_a})" if compound_a else name_a
            sub_from_b = f"({name_b})" if compound_b else name_b
            key_a, _, _ = winning_chain_with_key(full_graph, graph_a, {selenium_idx: _selanyl_prefix(sub_from_b)}, mol=mol)
            key_b, _, _ = winning_chain_with_key(full_graph, graph_b, {selenium_idx: _selanyl_prefix(sub_from_a)}, mol=mol)
            parent_root, sub_root = (n1, n2) if key_a <= key_b else (n2, n1)
    elif size1 > size2:
        parent_root, sub_root = n1, n2
    else:
        parent_root, sub_root = n2, n1

    sub_name, sub_compound = name_branch(full_graph, sub_root, selenium_idx, {}, mol=mol)
    if sub_compound:
        sub_name = f"({sub_name})"

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {selenium_idx: _selanyl_prefix(sub_name)}
    return name_from_carbon_graph(full_graph, parent_carbon_graph, terminals, mol=mol)
