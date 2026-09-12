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
branched R' substituent. Still out of scope: any unsaturation or non-
benzene ring, and any heteroatom other than the single sulfide sulfur (in
particular a disulfide S-S, or an oxidized sulfur -- sulfoxide/sulfone --
are separate functional groups, not in scope here).

- P-91.3/P-92: a molecule with
  one or more *specified* tetrahedral stereocenters on the parent (R)
  chain gets a "(<locant><R/S>,...)-" prefix, ascending locant order --
  same mechanism as `_ether.py`, using `_acyclic.py`'s
  `winning_chain_from_carbon_graph` for the parent chain's locant lookup.
  A stereocenter on the sulfanyl (R') substituent branch remains out of
  scope (raises `UnsupportedStructure`).

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_sulfide_chain`), mirroring `_ether.py`'s equivalent:
'sulfanyl' has no suffix form, so P-44.1.2.2 rule (1) makes the ring the
parent regardless of the other side's chain length. Both a direct
ring-sulfur bond and a chain spacer between the ring and the sulfide
sulfur are supported, confirmed via PubChem PUG REST for the underlying
structure/connectivity (CID 12144 'c1ccccc1SCC' -> 'ethylsulfanylbenzene',
'c1ccccc1CSCC' -> 'ethylsulfanylmethylbenzene', 'c1ccccc1SC(C)C' ->
'propan-2-ylsulfanylbenzene', 'c1ccccc1CSC(C)C' ->
'propan-2-ylsulfanylmethylbenzene') -- but, like `_ether.py`, a compound
R' side or a compound chain-spacer branch is still parenthesized here
even where PubChem's own auto-generated name omits the parentheses.
"""

from rdkit import Chem

from ._acyclic import longest_chain_length, winning_chain_from_carbon_graph, winning_chain_with_key
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    component_subgraph,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._substituents import name_branch


def _sulfanyl_prefix(name: str) -> str:
    return name + "sulfanyl"


def has_sulfide_shape(mol) -> bool:
    """True iff `mol` has exactly one sulfur atom, singly bonded to two
    carbons (a plain sulfide -S-, P-63.2.1) -- the shape this module
    accepts."""
    sulfurs = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if len(sulfurs) != 1:
        return False
    (sulfur,) = sulfurs
    return sulfur.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in sulfur.GetNeighbors())


def _validate_sulfide_atoms(mol, ring_atoms=frozenset()):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path (mirrors `_ether.py`'s
    `_validate_ether_atoms`): every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a non-aromatic chain carbon or the
    sulfide's own sulfur."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() not in (6, 16):
            raise UnsupportedStructure(
                "heteroatoms other than the sulfide sulfur are not "
                "supported yet (P-63.2.1 is restricted to a plain -S- "
                "sulfide here)"
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
            "alkenes/alkynes; not yet combined with a sulfide here)"
        )


def _name_benzene_ring_sulfide_chain(mol, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since the 'sulfanyl' prefix has no suffix form
    (module docstring), a single, otherwise-unsubstituted benzene ring is
    always the parent hydride, regardless of the other side's chain
    length. Handles both a direct ring-sulfur bond ('ethylsulfanylbenzene')
    and a chain spacer between the ring and the sulfide sulfur
    ('(2-ethylsulfanylethyl)benzene'), mirroring
    `_ether._name_benzene_ring_ether_chain` exactly (sulfur in place of
    oxygen, 'sulfanyl' in place of 'oxy')."""
    _validate_sulfide_atoms(mol, ring_atoms)
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter is not supported yet for a "
            "benzene-ring-substituent sulfide"
        )

    graph = adjacency(mol)
    (sulfur_idx,) = (idx for idx in graph if mol.GetAtomWithIdx(idx).GetAtomicNum() == 16)

    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a sulfide is not supported yet"
        )
    ring_atom, chain_root = attachment

    if chain_root == sulfur_idx:
        (r_prime,) = [n for n in graph[sulfur_idx] if n != ring_atom]
        sub_name, sub_compound = name_branch(graph, r_prime, sulfur_idx, {}, mol=mol)
        if sub_compound:
            sub_name = f"({sub_name})"
        return f"{_sulfanyl_prefix(sub_name)}benzene"

    blocked_graph = {node: [n for n in neighbors if n != sulfur_idx] for node, neighbors in graph.items()}
    del blocked_graph[sulfur_idx]
    reached, _ = bfs(blocked_graph, ring_atom)
    (r_prime,) = [n for n in graph[sulfur_idx] if n not in reached]

    sub_name, sub_compound = name_branch(graph, r_prime, sulfur_idx, {}, mol=mol)
    if sub_compound:
        sub_name = f"({sub_name})"
    sulfanyl_term = _sulfanyl_prefix(sub_name)
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, {sulfur_idx: sulfanyl_term}, mol=mol)
    if not is_compound:
        return f"{branch_name}benzene"
    if "(" in branch_name:
        return f"[{branch_name}]benzene"
    return f"({branch_name})benzene"


def name_sulfide(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_sulfide_chain(mol, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_sulfide_atoms(mol)

    (sulfur,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    sulfur_idx = sulfur.GetIdx()
    n1, n2 = (n.GetIdx() for n in sulfur.GetNeighbors())

    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    size1 = len(component_subgraph(carbon_graph, n1))
    size2 = len(component_subgraph(carbon_graph, n2))

    if size1 == size2:
        graph_a = component_subgraph(carbon_graph, n1)
        graph_b = component_subgraph(carbon_graph, n2)
        len_a, len_b = longest_chain_length(graph_a), longest_chain_length(graph_b)
        if len_a > len_b:
            parent_root, sub_root = n1, n2
        elif len_b > len_a:
            parent_root, sub_root = n2, n1
        else:
            name_a, compound_a = name_branch(full_graph, n1, sulfur_idx, {}, mol=mol)
            name_b, compound_b = name_branch(full_graph, n2, sulfur_idx, {}, mol=mol)
            sub_from_a = f"({name_a})" if compound_a else name_a
            sub_from_b = f"({name_b})" if compound_b else name_b
            key_a, _, _ = winning_chain_with_key(full_graph, graph_a, {sulfur_idx: _sulfanyl_prefix(sub_from_b)}, mol=mol)
            key_b, _, _ = winning_chain_with_key(full_graph, graph_b, {sulfur_idx: _sulfanyl_prefix(sub_from_a)}, mol=mol)
            parent_root, sub_root = (n1, n2) if key_a <= key_b else (n2, n1)
    elif size1 > size2:
        parent_root, sub_root = n1, n2
    else:
        parent_root, sub_root = n2, n1

    sub_name, sub_compound = name_branch(full_graph, sub_root, sulfur_idx, {}, mol=mol)
    if sub_compound:
        sub_name = f"({sub_name})"

    parent_carbon_graph = component_subgraph(carbon_graph, parent_root)
    terminals = {sulfur_idx: _sulfanyl_prefix(sub_name)}
    chain, name = winning_chain_from_carbon_graph(full_graph, parent_carbon_graph, terminals, mol=mol)

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
