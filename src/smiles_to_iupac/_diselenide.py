"""Naming of diselenides (the 'diselanyl' substituent prefix, R-Se-Se-R'),
restricted to two acyclic saturated hydrocarbon chains hung off a single
-Se-Se- bridge, per the IUPAC 2013 Recommendations ("the Blue Book"):

- Structurally the same shape as `_selenide.py`/`_sulfide.py`/
  `_peroxide.py`: one side (RH) is the parent hydride, the other
  (R'-Se-Se-) a substituent prefix on it, parent choice following the same
  P-44.3 skeletal-atom-count rule. Confirmed via PubChem PUG REST: CID
  23496 (`C[Se][Se]C`) -> "(methyldiselanyl)methane", CID 69405
  (`CC[Se][Se]CC`) -> "(ethyldiselanyl)ethane", CID 85591054
  (`CCC[Se][Se]C`) -> "1-(methyldiselanyl)propane" (a 3-carbon parent
  needs the locant, same rule as `_selenide.py`), CID 129678518
  (`C[Se][Se]CC`) -> "(methyldiselanyl)ethane" (2-carbon parent omits it).
- Unlike `_sulfide.py`/`_selenide.py`/`_peroxide.py`'s plain
  'alkylsulfanyl'/'alkylselanyl'/'alkylperoxy' prefixes, the R'-Se-Se-
  prefix here is *always* enclosed in parentheses, even when it isn't
  branched -- confirmed directly by every worked example above (contrast
  `_peroxide.py`'s parenthesis-free 'methylperoxymethane'). This appears
  to be P-16.3.3's disambiguation rule for a substituent name that itself
  begins with a syllable that looks like a multiplying prefix ('di-' in
  'diselanyl', which could otherwise be misread as a "di" multiplier on a
  plain 'selanyl' prefix). Because of this, this module can't reuse
  `_acyclic.py`'s `name_from_carbon_graph`/`terminals` shortcut the other
  three modules share (that shortcut always cites `is_compound=False`);
  it instead builds its own small parent-chain substituent list that
  forces the diselanyl prefix's `is_compound` flag to `True`, and its own
  locant-omission check (P-14.3.4.2(b) applies here regardless of the
  sole substituent being compound or not -- confirmed by the CID 129678518
  two-carbon case above, unlike `_acyclic.py`'s own equivalent check,
  which only omits the locant for a *non*-compound sole substituent; that
  narrower check is untouched here since no existing caller of
  `_acyclic.py` has ever exercised a compound substituent on a two-carbon
  chain).

- One side may also be a plain terminal -SeH instead of a second R'
  (R-Se-SeH): confirmed via PubChem PUG REST -- CID 101729611
  (`C[Se][SeH]`) -> "diselanylmethane" (mononuclear parent, no
  parentheses -- a *bare*, unprefixed 'diselanyl' name isn't itself a
  compound substituent, so P-16.3.3's disambiguation doesn't apply when
  no locant sits next to it), CID 173348765 (`CC[Se][SeH]`) ->
  "diselanylethane" (2-carbon parent, same no-locant/no-parens rule), CID
  174964680 (`CCC[Se][SeH]`) -> "1-(diselanyl)propane" (3-carbon parent --
  once a locant digit sits directly in front of 'diselanyl', the
  parentheses return, same as `_disulfide.py`'s identical case). Both
  seleniums terminal (H-Se-Se-H) stays out of scope -- no carbon parent to
  hang a name on.

Scope and out-of-scope structures are otherwise identical to
`_selenide.py`, selenium chained to selenium in place of a single
selenium -- see that module's docstring. In particular still out of
scope: a branched R' substituent, any unsaturation or ring, and any
heteroatom other than the diselenide's own two seleniums (a
monoselenide, or oxidized selenium, are separate functional groups).
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    group_substituents,
    is_plain_benzene_ring,
    longest_chains,
    lowest_locant_set,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_SELENIUM = 34
_BARE_TERMINAL_NAME = "diselanyl"


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_diselenide_shape(mol) -> bool:
    """True iff `mol` has exactly two selenium atoms, singly bonded to
    each other, where each selenium is additionally bonded either to
    exactly one carbon (R-Se-) or to nothing else (a terminal -SeH,
    filled by an implicit hydrogen) -- but not both terminal at once
    (H-Se-Se-H has no carbon parent to hang a name on, out of scope)."""
    seleniums = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _SELENIUM]
    if len(seleniums) != 2:
        return False
    se1, se2 = seleniums
    bond = mol.GetBondBetweenAtoms(se1.GetIdx(), se2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return False
    if se1.GetDegree() not in (1, 2) or se2.GetDegree() not in (1, 2):
        return False
    if se1.GetDegree() == 1 and se1.GetTotalNumHs() != 1:
        return False
    if se2.GetDegree() == 1 and se2.GetTotalNumHs() != 1:
        return False
    if se1.GetDegree() == 1 and se2.GetDegree() == 1:
        return False
    others = []
    for se, other_se in ((se1, se2), (se2, se1)):
        if se.GetDegree() == 2:
            (other,) = [n for n in se.GetNeighbors() if n.GetIdx() != other_se.GetIdx()]
            others.append(other)
    return all(o.GetAtomicNum() == 6 for o in others)


def _validate_and_find_diselenide(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom/ring/unsaturation rejections below so
    `name_diselenide`'s benzene-ring-substituent path (see
    `_name_benzene_ring_diselenide_chain`) can reuse this same validation
    for the rest of the molecule. Empty by default, so every other
    caller's behavior is unchanged."""
    if not has_diselenide_shape(mol):
        raise UnsupportedStructure(
            "no plain diselenide (R-Se-Se-R') skeleton found; this module "
            "only handles diselenides"
        )
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, _SELENIUM):
            raise UnsupportedStructure(
                "heteroatoms other than the diselenide's own two seleniums "
                "are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if any(a not in aromatic_ring_atoms or b not in aromatic_ring_atoms for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see P-31 for "
            "alkenes/alkynes; not yet combined with a diselenide here)"
        )
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0 and not (
        ring_info.NumRings() == 1 and set(ring_info.AtomRings()[0]) == set(aromatic_ring_atoms)
    ):
        raise UnsupportedStructure(
            "rings are not supported yet, other than the separate "
            "benzene-ring-substituent path"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    se1, se2 = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _SELENIUM)
    c1 = None
    if se1.GetDegree() == 2:
        (c1,) = [n.GetIdx() for n in se1.GetNeighbors() if n.GetIdx() != se2.GetIdx()]
    c2 = None
    if se2.GetDegree() == 2:
        (c2,) = [n.GetIdx() for n in se2.GetNeighbors() if n.GetIdx() != se1.GetIdx()]
    return se1.GetIdx(), se2.GetIdx(), c1, c2


def _substituents_for_chain(graph, chain, terminals, mol=None):
    """Like `_acyclic._substituents_for_chain`, but a `terminals` leaf is
    always cited as a compound (parenthesized) substituent -- see module
    docstring for why the diselanyl prefix needs this, unlike the plain
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
                entries.append(name_branch(graph, root, atom, {}, mol=mol))
        substituents[position] = entries
    return substituents


def _name_from_substituents(chain_length, grouped):
    if chain_length == 1 and grouped:
        (name,) = grouped
        if name == _BARE_TERMINAL_NAME:
            # No locant sits next to it here, so the bare (unprefixed)
            # name needs no P-16.3.3 parentheses -- confirmed by
            # 'diselanylmethane' (CID 101729611), unlike a genuinely
            # compound name like 'methyldiselanyl' (still parenthesized
            # below via the ordinary path).
            return name + alkane_name(chain_length)
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(chain_length)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if chain_length == 2 and total_subs == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly
        # one substituent has only one possible structure regardless of
        # numbering direction, so the locant is omittable -- this holds
        # whether or not that substituent is itself compound (confirmed
        # by 'diselanyl', unlike `_acyclic.py`'s own narrower check). A
        # bare (unprefixed) 'diselanyl' name is the exception -- no locant
        # here either, so no parentheses (CID 173348765
        # 'diselanylethane'), unlike the >=3-carbon case below where a
        # locant digit does sit next to it (CID 174964680
        # '1-(diselanyl)propane').
        (name,) = grouped
        if name == _BARE_TERMINAL_NAME:
            display_name = name
        else:
            display_name = f"({name})" if grouped[name]["compound"] else name
        return display_name + alkane_name(chain_length)
    prefix = format_substituent_prefixes(grouped)
    return prefix + alkane_name(chain_length)


def _candidate_key(chain_length, substituents):
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, grouped)
    return (-total_count, locant_set, citation_locants, name), name


def _name_parent_chain(full_graph, carbon_graph, terminals, mol=None):
    chains = longest_chains(carbon_graph)
    chain_length = len(chains[0])

    best_key = None
    best_name = None
    best_chain = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            substituents = _substituents_for_chain(full_graph, candidate, terminals, mol=mol)
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_chain = key, name, candidate
    return best_chain, best_name


def _name_benzene_ring_diselenide_chain(mol, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since 'diselanyl' has no suffix form (module
    docstring), a single, otherwise-unsubstituted benzene ring is senior
    to a chain of the same (plain-hydrocarbon) class regardless of the
    chain's length -- mirrors `_disulfide.py`'s
    `_name_benzene_ring_disulfide_chain`. The ring is always the parent
    hydride; the whole R-Se-Se- fragment hung *directly* off the ring is
    a single '...diselanyl' substituent prefix on it (a chain spacer
    between the ring and the near selenium is out of scope, same reason
    as `_disulfide.py`). Confirmed via PubChem: CID 59041517
    ('(methyldiselanyl)benzene'), CID 71327844
    ('(ethyldiselanyl)benzene')."""
    se1_idx, se2_idx, c1, c2 = _validate_and_find_diselenide(mol, aromatic_ring_atoms=ring_atoms)
    if c1 is None or c2 is None:
        raise UnsupportedStructure(
            "a -SeH terminal combined with a benzene-ring-substituent "
            "diselenide is out of scope for this module"
        )
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter is not supported yet for a "
            "benzene-ring-substituent diselenide"
        )

    full_graph = adjacency(mol)
    attachment = ring_chain_attachment(full_graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a diselenide chain is not supported yet"
        )
    ring_atom, root = attachment

    if root == se1_idx:
        other_se, other_root = se2_idx, c2
    elif root == se2_idx:
        other_se, other_root = se1_idx, c1
    else:
        raise UnsupportedStructure(
            "a chain spacer between the benzene ring and the diselenide's "
            "near selenium is not supported yet (only a direct ring-"
            "selenium bond is, see module docstring)"
        )

    sub_name, sub_compound = name_branch(full_graph, other_root, other_se, {}, mol=mol)
    if sub_compound:
        raise UnsupportedStructure(
            "a branched alkyldiselanyl substituent is not supported yet"
        )
    return f"({sub_name}diselanyl)benzene"


def name_diselenide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_diselenide_chain(mol, ring_atoms)
    se1_idx, se2_idx, c1, c2 = _validate_and_find_diselenide(mol)
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)

    if c1 is None or c2 is None:
        parent_root, parent_se = (c2, se2_idx) if c1 is None else (c1, se1_idx)
        terminals = {parent_se: _BARE_TERMINAL_NAME}
    else:
        size1 = len(_component_subgraph(carbon_graph, c1))
        size2 = len(_component_subgraph(carbon_graph, c2))

        if size1 == size2:
            _, compound_a = name_branch(full_graph, c1, se1_idx, {}, mol=mol)
            _, compound_b = name_branch(full_graph, c2, se2_idx, {}, mol=mol)
            if compound_a and compound_b:
                raise UnsupportedStructure(
                    "a diselenide tied in skeletal-atom count with both sides "
                    "branched is not supported yet"
                )
            if compound_a:
                parent_root, parent_se, sub_root, sub_se = c2, se2_idx, c1, se1_idx
            else:
                parent_root, parent_se, sub_root, sub_se = c1, se1_idx, c2, se2_idx
        elif size1 > size2:
            parent_root, parent_se, sub_root, sub_se = c1, se1_idx, c2, se2_idx
        else:
            parent_root, parent_se, sub_root, sub_se = c2, se2_idx, c1, se1_idx

        sub_name, sub_compound = name_branch(full_graph, sub_root, sub_se, {}, mol=mol)
        if sub_compound:
            raise UnsupportedStructure(
                "a branched alkyldiselanyl substituent is not supported yet"
            )
        terminals = {parent_se: sub_name + "diselanyl"}

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    chain, name = _name_parent_chain(full_graph, parent_carbon_graph, terminals, mol=mol)

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter must lie on the winning parent
    # chain; one on the diselanyl (R') substituent branch is out of
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
