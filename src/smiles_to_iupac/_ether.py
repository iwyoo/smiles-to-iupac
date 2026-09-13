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
- P-91.3/P-92: a molecule with one
  or more *specified* tetrahedral stereocenters on the parent (R) chain
  gets a "(<locant><R/S>,...)-" prefix, ascending locant order, e.g.
  '(2S)-2-ethoxybutane' -- same mechanism as `_acetal.py`, using
  `_acyclic.py`'s `winning_chain_from_carbon_graph` for the parent
  chain's locant lookup. A stereocenter on the alkoxy (R') substituent
  branch remains out of scope (raises `UnsupportedStructure`).
- A tie in total skeletal-atom count between the two chains is resolved in
  two steps, both reusing `_acyclic.py`'s existing P-44.3/P-45.2 machinery
  rather than a new algorithm: first by each side's own longest achievable
  chain length (`longest_chain_length` -- e.g. a neopentyl arm's 5 carbons
  max out at chain length 3 around its quaternary carbon, while an
  isopentyl arm's 5 carbons reach chain length 4, so isopentyl wins
  outright, PubChem CID-confirmed 'CC(C)(C)COCCC(C)C' ->
  '1-(2,2-dimethylpropoxy)-3-methylbutane'); only if that also ties (both
  sides reduce to the same chain length, e.g. an isobutyl arm vs. a
  tert-butyl arm, both length 3) does the choice fall to comparing each
  side's resulting locant set as parent (`winning_chain_with_key`,
  PubChem CID-confirmed 'CC(C)COC(C)(C)C' ->
  '2-methyl-1-[(2-methylpropan-2-yl)oxy]propane': isobutyl-as-parent gets
  locants {1,2}, tert-butyl-as-parent gets {2,2}, so isobutyl wins).

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one oxygen, or an oxygen not shaped like a plain ether (degree
  != 2, a non-carbon neighbor) — routed to a different module by `core.py`
  before this one is even tried.
- Any unsaturation, any ring other than a single benzene ring (see below),
  or any heteroatom other than the single ether oxygen.

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_ether_chain`), mirroring `_nitro.py`'s equivalent:
'oxy' has no suffix form, so P-44.1.2.2 rule (1) makes the ring the parent
regardless of the other side's chain length. Both a direct ring-oxygen
bond (P-63.2.2.1.1's own 'alkoxybenzene' shape) and a chain spacer between
the ring and the ether oxygen are supported, confirmed via PubChem PUG
REST for the underlying structure/connectivity (CID 7500 'c1ccccc1OCC' ->
'ethoxybenzene', 'c1ccccc1COCC' -> 'ethoxymethylbenzene', 'c1ccccc1OC(C)C'
-> 'propan-2-yloxybenzene', 'c1ccccc1COC(C)C' ->
'propan-2-yloxymethylbenzene') — but, like `_nitro.py`/`_azide.py`
(see their own PubChem-vs-PIN discrepancy notes in `test_nitro.py`/
`test_azide.py`), a compound R' side or a compound chain-spacer branch is
still parenthesized here even where PubChem's own auto-generated name
omits the parentheses, escalating to square brackets rather than nesting
round ones when a '(...)oxy' term already sits inside the branch.
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

_CONTRACTED_OXY = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy"}


def _oxy_prefix(name: str) -> str:
    return _CONTRACTED_OXY.get(name, name + "oxy")


def has_ether_shape(mol) -> bool:
    """True iff `mol` has exactly one oxygen atom, singly bonded to two
    carbons (a plain ether -O-, P-63.2.1) — the shape this module accepts.
    A furan ring's own aromatic oxygen is excluded from the count so it
    doesn't misclaim a molecule whose real functional group lies
    elsewhere (P-616 M4)."""
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and not atom.GetIsAromatic()]
    if len(oxygens) != 1:
        return False
    (oxygen,) = oxygens
    return oxygen.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in oxygen.GetNeighbors())


def _validate_ether_atoms(mol, ring_atoms=frozenset()):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path (mirrors `_nitro.py`'s
    `_validate_nitro_atoms`): every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a non-aromatic chain carbon or the
    ether's own oxygen."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() not in (6, 8):
            raise UnsupportedStructure(
                "heteroatoms other than the ether oxygen are not supported "
                "yet (P-63.2.1 is restricted to a plain -O- ether here)"
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
            "alkenes/alkynes; not yet combined with an ether here)"
        )


def _name_benzene_ring_ether_chain(mol, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since the 'oxy' prefix has no suffix form
    (module docstring), a single, otherwise-unsubstituted benzene ring is
    always the parent hydride, regardless of the other side's chain
    length. Handles both a direct ring-oxygen bond ('ethoxybenzene') and
    a chain spacer between the ring and the ether oxygen
    ('(2-ethoxyethyl)benzene').

    A branched R' still gets the plain two-chain path's own parenthesized
    '(...)oxy' treatment (P-63.2.2.1.1's worked example, 'butan-2-yl' ->
    '(butan-2-yl)oxy'). A chain-spacer branch also follows this project's
    existing benzene-ring-chain convention (`_nitro.py`/`_azide.py`): the
    whole branch is parenthesized whenever it's compound, even though
    PubChem's own auto-generated name sometimes omits those parentheses
    (confirmed for azide, see `test_azide.py`) -- escalating to square
    brackets when the branch name already contains a '(...)oxy' term of
    its own, to avoid nesting round brackets (mirrors the plain two-chain
    path's own bracket-escalation example in the module docstring)."""
    _validate_ether_atoms(mol, ring_atoms)
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter is not supported yet for a "
            "benzene-ring-substituent ether"
        )

    graph = adjacency(mol)
    (oxygen_idx,) = (idx for idx in graph if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)

    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside an ether is not supported yet"
        )
    ring_atom, chain_root = attachment

    if chain_root == oxygen_idx:
        (r_prime,) = [n for n in graph[oxygen_idx] if n != ring_atom]
        sub_name, sub_compound = name_branch(graph, r_prime, oxygen_idx, {}, mol=mol)
        oxy_term = _oxy_prefix(sub_name)
        if sub_compound:
            oxy_term = f"({oxy_term})"
        return f"{oxy_term}benzene"

    blocked_graph = {node: [n for n in neighbors if n != oxygen_idx] for node, neighbors in graph.items()}
    del blocked_graph[oxygen_idx]
    reached, _ = bfs(blocked_graph, ring_atom)
    (r_prime,) = [n for n in graph[oxygen_idx] if n not in reached]

    sub_name, sub_compound = name_branch(graph, r_prime, oxygen_idx, {}, mol=mol)
    oxy_term = _oxy_prefix(sub_name)
    if sub_compound:
        oxy_term = f"({oxy_term})"
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, {oxygen_idx: oxy_term}, mol=mol)
    if not is_compound:
        return f"{branch_name}benzene"
    if "(" in branch_name:
        return f"[{branch_name}]benzene"
    return f"({branch_name})benzene"


def name_ether(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_ether_chain(mol, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_ether_atoms(mol)

    (oxygen,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    oxygen_idx = oxygen.GetIdx()
    n1, n2 = (n.GetIdx() for n in oxygen.GetNeighbors())

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
            name_a, compound_a = name_branch(full_graph, n1, oxygen_idx, {}, mol=mol)
            name_b, compound_b = name_branch(full_graph, n2, oxygen_idx, {}, mol=mol)
            oxy_from_a = _oxy_prefix(name_a)
            oxy_from_b = _oxy_prefix(name_b)
            if compound_a:
                oxy_from_a = f"({oxy_from_a})"
            if compound_b:
                oxy_from_b = f"({oxy_from_b})"
            key_a, _, _ = winning_chain_with_key(full_graph, graph_a, {oxygen_idx: oxy_from_b}, mol=mol)
            key_b, _, _ = winning_chain_with_key(full_graph, graph_b, {oxygen_idx: oxy_from_a}, mol=mol)
            parent_root, sub_root = (n1, n2) if key_a <= key_b else (n2, n1)
    elif size1 > size2:
        parent_root, sub_root = n1, n2
    else:
        parent_root, sub_root = n2, n1

    sub_name, sub_compound = name_branch(full_graph, sub_root, oxygen_idx, {}, mol=mol)
    oxy_term = _oxy_prefix(sub_name)
    if sub_compound:
        oxy_term = f"({oxy_term})"

    parent_carbon_graph = component_subgraph(carbon_graph, parent_root)
    terminals = {oxygen_idx: oxy_term}
    chain, name = winning_chain_from_carbon_graph(full_graph, parent_carbon_graph, terminals, mol=mol)

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
