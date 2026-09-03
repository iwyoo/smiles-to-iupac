"""Naming of disulfides (the 'disulfanyl' substituent prefix, R-S-S-R'),
restricted to two acyclic saturated hydrocarbon chains hung off a single
-S-S- bridge, per the IUPAC 2013 Recommendations ("the Blue Book"):

- Structurally the same shape as `_diselenide.py`/`_ditelluride.py` (their
  own docstrings cover the shared reasoning in full): one side (RH) is
  the parent hydride, the other (R'-S-S-) a substituent prefix on it,
  always enclosed in parentheses (P-16.3.3's disambiguation rule for a
  substituent name beginning with a syllable that looks like a
  multiplying prefix, 'di-' in 'disulfanyl'), which is why this module
  can't reuse `_acyclic.py`'s `name_from_carbon_graph`/`terminals`
  shortcut `_sulfide.py` itself uses and instead mirrors
  `_diselenide.py`'s own bespoke parent-chain substituent/locant logic,
  sulfur in place of selenium. `_thiol.py`/`_sulfide.py` both explicitly
  scoped a disulfide -S-S- out as future work; this module is that
  follow-up. Confirmed via PubChem PUG REST -- all four cases (unlike
  `_diselenide.py`/`_ditelluride.py`, every one of these has a real
  PubChem-listed compound): CID 12232 (`CSSC`) ->
  "(methyldisulfanyl)methane", CID 8077 (`CCSSCC`) ->
  "(ethyldisulfanyl)ethane", CID 123388 (`CSSCC`) ->
  "(methyldisulfanyl)ethane" (2-carbon parent omits the locant, same
  P-14.3.4.2(b) rule as `_diselenide.py`'s identical case), CID 16592
  (`CCCSSC`) -> "1-(methyldisulfanyl)propane" (3-carbon parent needs it).

- One side may also be a plain terminal -SH instead of a second R'
  (a "perthiol"/hydrodisulfide, R-S-SH): confirmed via PubChem PUG REST --
  CID 522059 (`CSS`) -> "disulfanylmethane" (mononuclear parent, no
  parentheses -- unlike the R-S-S-R' case above, a *bare*, unprefixed
  'disulfanyl' name isn't itself a compound substituent, so P-16.3.3's
  disambiguation doesn't apply when no locant sits next to it), CID 94671
  (`CCSS`) -> "disulfanylethane" (2-carbon parent, same no-locant/no-parens
  rule), CID 6428842 (`CCCSS`) -> "1-(disulfanyl)propane" (3-carbon parent
  -- once a locant digit sits directly in front of 'disulfanyl', the
  parentheses return to keep the digit from misreading as part of the
  'di-' syllable). Both sulfurs terminal (H-S-S-H, disulfane itself) stays
  out of scope -- no carbon parent to hang a name on.

Scope and out-of-scope structures are otherwise identical to
`_diselenide.py`, sulfur chained to sulfur in place of selenium chained
to selenium. In particular still out of scope: a branched R'
substituent, any unsaturation or ring, and any heteroatom other than the
disulfide's own two sulfurs (a plain sulfide, thiol, or oxidized sulfur,
are separate functional groups).
"""

from rdkit import Chem

from ._acyclic import _longest_chains
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    lowest_locant_set,
    non_single_bonds,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_SULFUR = 16
_BARE_TERMINAL_NAME = "disulfanyl"


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_disulfide_shape(mol) -> bool:
    """True iff `mol` has exactly two sulfur atoms, singly bonded to each
    other, where each sulfur is additionally bonded either to exactly one
    carbon (R-S-) or to nothing else (a terminal -SH, filled by an
    implicit hydrogen) -- but not both terminal at once (H-S-S-H has no
    carbon parent to hang a name on, out of scope)."""
    sulfurs = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _SULFUR]
    if len(sulfurs) != 2:
        return False
    s1, s2 = sulfurs
    bond = mol.GetBondBetweenAtoms(s1.GetIdx(), s2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return False
    if s1.GetDegree() not in (1, 2) or s2.GetDegree() not in (1, 2):
        return False
    if s1.GetDegree() == 1 and s1.GetTotalNumHs() != 1:
        return False
    if s2.GetDegree() == 1 and s2.GetTotalNumHs() != 1:
        return False
    if s1.GetDegree() == 1 and s2.GetDegree() == 1:
        return False
    others = []
    for s, other_s in ((s1, s2), (s2, s1)):
        if s.GetDegree() == 2:
            (other,) = [n for n in s.GetNeighbors() if n.GetIdx() != other_s.GetIdx()]
            others.append(other)
    return all(o.GetAtomicNum() == 6 for o in others)


def _validate_and_find_disulfide(mol):
    if not has_disulfide_shape(mol):
        raise UnsupportedStructure(
            "no plain disulfide (R-S-S-R') skeleton found; this module "
            "only handles disulfides"
        )
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, _SULFUR):
            raise UnsupportedStructure(
                "heteroatoms other than the disulfide's own two sulfurs "
                "are not supported yet"
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
            "alkenes/alkynes; not yet combined with a disulfide here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    s1, s2 = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _SULFUR)
    c1 = None
    if s1.GetDegree() == 2:
        (c1,) = [n.GetIdx() for n in s1.GetNeighbors() if n.GetIdx() != s2.GetIdx()]
    c2 = None
    if s2.GetDegree() == 2:
        (c2,) = [n.GetIdx() for n in s2.GetNeighbors() if n.GetIdx() != s1.GetIdx()]
    return s1.GetIdx(), s2.GetIdx(), c1, c2


def _substituents_for_chain(graph, chain, terminals):
    """Like `_acyclic._substituents_for_chain`, but a `terminals` leaf is
    always cited as a compound (parenthesized) substituent -- see module
    docstring for why the disulfanyl prefix needs this, unlike the plain
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
        (name,) = grouped
        if name == _BARE_TERMINAL_NAME:
            # No locant sits next to it here, so the bare (unprefixed)
            # name needs no P-16.3.3 parentheses -- confirmed by
            # 'disulfanylmethane' (CID 522059), unlike a genuinely
            # compound name like 'methyldisulfanyl' (still parenthesized
            # below via the ordinary path).
            return name + alkane_name(chain_length)
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(chain_length)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if chain_length == 2 and total_subs == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly
        # one substituent has only one possible structure regardless of
        # numbering direction, so the locant is omittable -- this holds
        # whether or not that substituent is itself compound (confirmed
        # by 'disulfanyl', same as `_diselenide.py`'s identical rule).
        # A bare (unprefixed) 'disulfanyl' name is the exception -- no
        # locant here either, so no parentheses (CID 94671
        # 'disulfanylethane'), unlike the >=3-carbon case below where a
        # locant digit does sit next to it (CID 6428842
        # '1-(disulfanyl)propane').
        (name,) = grouped
        if name == _BARE_TERMINAL_NAME:
            display_name = name
        else:
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
    best_chain = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            substituents = _substituents_for_chain(full_graph, candidate, terminals)
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_chain = key, name, candidate
    return best_chain, best_name


def name_disulfide(mol) -> str:
    s1_idx, s2_idx, c1, c2 = _validate_and_find_disulfide(mol)
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)

    if c1 is None or c2 is None:
        parent_root, parent_s = (c2, s2_idx) if c1 is None else (c1, s1_idx)
        terminals = {parent_s: _BARE_TERMINAL_NAME}
    else:
        size1 = len(_component_subgraph(carbon_graph, c1))
        size2 = len(_component_subgraph(carbon_graph, c2))

        if size1 == size2:
            _, compound_a = name_branch(full_graph, c1, s1_idx, {})
            _, compound_b = name_branch(full_graph, c2, s2_idx, {})
            if compound_a and compound_b:
                raise UnsupportedStructure(
                    "a disulfide tied in skeletal-atom count with both sides "
                    "branched is not supported yet"
                )
            if compound_a:
                parent_root, parent_s, sub_root, sub_s = c2, s2_idx, c1, s1_idx
            else:
                parent_root, parent_s, sub_root, sub_s = c1, s1_idx, c2, s2_idx
        elif size1 > size2:
            parent_root, parent_s, sub_root, sub_s = c1, s1_idx, c2, s2_idx
        else:
            parent_root, parent_s, sub_root, sub_s = c2, s2_idx, c1, s1_idx

        sub_name, sub_compound = name_branch(full_graph, sub_root, sub_s, {})
        if sub_compound:
            raise UnsupportedStructure(
                "a branched alkyldisulfanyl substituent is not supported yet"
            )
        terminals = {parent_s: sub_name + "disulfanyl"}

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    chain, name = _name_parent_chain(full_graph, parent_carbon_graph, terminals)

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter must lie on the winning parent
    # chain; one on the disulfanyl (R') substituent branch is out of
    # scope, mirroring `_sulfide.py`'s identical treatment.
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    if any(atom not in position_of for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the "
            "principal chain is not supported yet (see P-92)"
        )
    labels = sorted((position_of[atom], code) for atom, code in stereo)
    prefix = ",".join(f"{locant}{code}" for locant, code in labels)
    return f"({prefix})-{name}"
