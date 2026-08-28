"""Naming of ditellurides (the 'ditellanyl' substituent prefix,
R-Te-Te-R'), restricted to two acyclic saturated hydrocarbon chains hung
off a single -Te-Te- bridge, per the IUPAC 2013 Recommendations ("the
Blue Book"):

- Structurally the same shape as `_diselenide.py` (its own docstring
  covers the shared reasoning in full): one side (RH) is the parent
  hydride, the other (R'-Te-Te-) a substituent prefix on it, always
  enclosed in parentheses (P-16.3.3's disambiguation rule for a
  substituent name beginning with a syllable that looks like a
  multiplying prefix, 'di-' in 'ditellanyl'), which is why this module
  can't reuse `_acyclic.py`'s `name_from_carbon_graph`/`terminals`
  shortcut either and instead mirrors `_diselenide.py`'s own bespoke
  parent-chain substituent/locant logic, tellurium in place of selenium.
  Confirmed via PubChem PUG REST: CID 88493 (`C[Te][Te]C`) ->
  "(methylditellanyl)methane", CID 141264 (`CC[Te][Te]CC`) ->
  "(ethylditellanyl)ethane", CID 86011476 (`C[Te][Te]CC`) ->
  "(methylditellanyl)ethane" (2-carbon parent omits the locant, same
  P-14.3.4.2(b) rule as `_diselenide.py`'s identical case). The analogous
  3-carbon-parent locant case (`CCC[Te][Te]C`) has no PubChem-listed
  compound to independently confirm (CID 0), so it's included as a
  reviewed (eyeballed), not independently verified, single-axis extension
  of the already-confirmed `_diselenide.py` rule for that exact shape.

Scope and out-of-scope structures are otherwise identical to
`_diselenide.py`, tellurium chained to tellurium in place of selenium
chained to selenium. In particular still out of scope: a branched R'
substituent, any unsaturation or ring, and any heteroatom other than the
ditelluride's own two telluriums (a monotelluride, or oxidized tellurium,
are separate functional groups).
"""

from rdkit import Chem

from ._acyclic import _longest_chains
from ._common import UnsupportedStructure, adjacency, bfs, carbon_adjacency, lowest_locant_set, non_single_bonds
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_TELLURIUM = 52


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_ditelluride_shape(mol) -> bool:
    """True iff `mol` has exactly two tellurium atoms, bonded to each
    other (single bond) and each also singly bonded to exactly one carbon
    (a plain R-Te-Te-R' ditelluride) -- the shape this module accepts."""
    telluriums = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _TELLURIUM]
    if len(telluriums) != 2:
        return False
    te1, te2 = telluriums
    bond = mol.GetBondBetweenAtoms(te1.GetIdx(), te2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return False
    if te1.GetDegree() != 2 or te2.GetDegree() != 2:
        return False
    (other1,) = [n for n in te1.GetNeighbors() if n.GetIdx() != te2.GetIdx()]
    (other2,) = [n for n in te2.GetNeighbors() if n.GetIdx() != te1.GetIdx()]
    return other1.GetAtomicNum() == 6 and other2.GetAtomicNum() == 6


def _validate_and_find_ditelluride(mol):
    if not has_ditelluride_shape(mol):
        raise UnsupportedStructure(
            "no plain ditelluride (R-Te-Te-R') skeleton found; this module "
            "only handles ditellurides"
        )
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, _TELLURIUM):
            raise UnsupportedStructure(
                "heteroatoms other than the ditelluride's own two "
                "telluriums are not supported yet"
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
            "alkenes/alkynes; not yet combined with a ditelluride here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    te1, te2 = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _TELLURIUM)
    (c1,) = [n.GetIdx() for n in te1.GetNeighbors() if n.GetIdx() != te2.GetIdx()]
    (c2,) = [n.GetIdx() for n in te2.GetNeighbors() if n.GetIdx() != te1.GetIdx()]
    return te1.GetIdx(), te2.GetIdx(), c1, c2


def _substituents_for_chain(graph, chain, terminals):
    """Like `_acyclic._substituents_for_chain`, but a `terminals` leaf is
    always cited as a compound (parenthesized) substituent -- see module
    docstring for why the ditellanyl prefix needs this, unlike the plain
    `terminals` shortcut `_sulfide.py`/`_selenide.py`/`_peroxide.py`
    share via `_acyclic.name_from_carbon_graph`."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set]
        if not branch_roots:
            continue
        entries = []
        for root in branch_roots:
            if root in terminals:
                entries.append((terminals[root], True))
            else:
                entries.append(name_branch(graph, root, atom, {}))
        substituents[position] = entries
    return substituents


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, grouped):
    if chain_length == 1 and grouped:
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(chain_length)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if chain_length == 2 and total_subs == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly
        # one substituent has only one possible structure regardless of
        # numbering direction, so the locant is omittable -- this holds
        # whether or not that substituent is itself compound (confirmed
        # by 'ditellanyl', same as `_diselenide.py`'s identical rule).
        (name,) = grouped
        display_name = f"({name})" if grouped[name]["compound"] else name
        return display_name + alkane_name(chain_length)
    prefix = format_substituent_prefixes(grouped)
    return prefix + alkane_name(chain_length)


def _candidate_key(chain_length, substituents):
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, grouped)
    return (-total_count, locant_set, citation_locants, name), name


def _name_parent_chain(full_graph, carbon_graph, terminals):
    chains = _longest_chains(carbon_graph)
    chain_length = len(chains[0])

    best_key = None
    best_name = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            substituents = _substituents_for_chain(full_graph, candidate, terminals)
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_ditelluride(mol) -> str:
    te1_idx, te2_idx, c1, c2 = _validate_and_find_ditelluride(mol)
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    size1 = len(_component_subgraph(carbon_graph, c1))
    size2 = len(_component_subgraph(carbon_graph, c2))

    if size1 == size2:
        _, compound_a = name_branch(full_graph, c1, te1_idx, {})
        _, compound_b = name_branch(full_graph, c2, te2_idx, {})
        if compound_a and compound_b:
            raise UnsupportedStructure(
                "a ditelluride tied in skeletal-atom count with both sides "
                "branched is not supported yet"
            )
        if compound_a:
            parent_root, parent_te, sub_root, sub_te = c2, te2_idx, c1, te1_idx
        else:
            parent_root, parent_te, sub_root, sub_te = c1, te1_idx, c2, te2_idx
    elif size1 > size2:
        parent_root, parent_te, sub_root, sub_te = c1, te1_idx, c2, te2_idx
    else:
        parent_root, parent_te, sub_root, sub_te = c2, te2_idx, c1, te1_idx

    sub_name, sub_compound = name_branch(full_graph, sub_root, sub_te, {})
    if sub_compound:
        raise UnsupportedStructure(
            "a branched alkylditellanyl substituent is not supported yet"
        )

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    terminals = {parent_te: sub_name + "ditellanyl"}
    return _name_parent_chain(full_graph, parent_carbon_graph, terminals)
