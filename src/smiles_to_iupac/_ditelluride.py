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

- One side may also be a plain terminal -TeH instead of a second R'
  (R-Te-TeH): confirmed via PubChem PUG REST -- CID 173350364
  (`CC[Te][TeH]`) -> "ditellanylethane" (2-carbon parent, no
  parentheses -- a *bare*, unprefixed 'ditellanyl' name isn't itself a
  compound substituent, so P-16.3.3's disambiguation doesn't apply when
  no locant sits next to it, mirroring `_disulfide.py`'s/
  `_diselenide.py`'s identical rule). The analogous 1-carbon and
  3-carbon-parent cases (`C[Te][TeH]`, `CCC[Te][TeH]`) have no
  PubChem-listed compound to independently confirm (CID 0), so they're
  included as a reviewed (eyeballed), not independently verified,
  single-axis extension of the already-confirmed `_disulfide.py`/
  `_diselenide.py` rule for that exact shape. Both telluriums terminal
  (H-Te-Te-H) stays out of scope -- no carbon parent to hang a name on.

Scope and out-of-scope structures are otherwise identical to
`_diselenide.py`, tellurium chained to tellurium in place of selenium
chained to selenium. In particular still out of scope: a branched R'
substituent, any unsaturation or ring, and any heteroatom other than the
ditelluride's own two telluriums (a monotelluride, or oxidized tellurium,
are separate functional groups).
"""

from rdkit import Chem

from ._acyclic import _longest_chains
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    is_plain_benzene_ring,
    lowest_locant_set,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_TELLURIUM = 52
_BARE_TERMINAL_NAME = "ditellanyl"


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def has_ditelluride_shape(mol) -> bool:
    """True iff `mol` has exactly two tellurium atoms, singly bonded to
    each other, where each tellurium is additionally bonded either to
    exactly one carbon (R-Te-) or to nothing else (a terminal -TeH,
    filled by an implicit hydrogen) -- but not both terminal at once
    (H-Te-Te-H has no carbon parent to hang a name on, out of scope)."""
    telluriums = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _TELLURIUM]
    if len(telluriums) != 2:
        return False
    te1, te2 = telluriums
    bond = mol.GetBondBetweenAtoms(te1.GetIdx(), te2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return False
    if te1.GetDegree() not in (1, 2) or te2.GetDegree() not in (1, 2):
        return False
    if te1.GetDegree() == 1 and te1.GetTotalNumHs() != 1:
        return False
    if te2.GetDegree() == 1 and te2.GetTotalNumHs() != 1:
        return False
    if te1.GetDegree() == 1 and te2.GetDegree() == 1:
        return False
    others = []
    for te, other_te in ((te1, te2), (te2, te1)):
        if te.GetDegree() == 2:
            (other,) = [n for n in te.GetNeighbors() if n.GetIdx() != other_te.GetIdx()]
            others.append(other)
    return all(o.GetAtomicNum() == 6 for o in others)


def _validate_and_find_ditelluride(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom/ring/unsaturation rejections below so
    `name_ditelluride`'s benzene-ring-substituent path (see
    `_name_benzene_ring_ditelluride_chain`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged."""
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
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if any(a not in aromatic_ring_atoms or b not in aromatic_ring_atoms for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see P-31 for "
            "alkenes/alkynes; not yet combined with a ditelluride here)"
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

    te1, te2 = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _TELLURIUM)
    c1 = None
    if te1.GetDegree() == 2:
        (c1,) = [n.GetIdx() for n in te1.GetNeighbors() if n.GetIdx() != te2.GetIdx()]
    c2 = None
    if te2.GetDegree() == 2:
        (c2,) = [n.GetIdx() for n in te2.GetNeighbors() if n.GetIdx() != te1.GetIdx()]
    return te1.GetIdx(), te2.GetIdx(), c1, c2


def _substituents_for_chain(graph, chain, terminals, mol=None):
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
                entries.append(name_branch(graph, root, atom, {}, mol=mol))
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
            # name needs no P-16.3.3 parentheses -- mirrors
            # `_disulfide.py`'s/`_diselenide.py`'s confirmed identical
            # rule ('methylditellanyl' etc. is still parenthesized below
            # via the ordinary path).
            return name + alkane_name(chain_length)
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(chain_length)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if chain_length == 2 and total_subs == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly
        # one substituent has only one possible structure regardless of
        # numbering direction, so the locant is omittable -- this holds
        # whether or not that substituent is itself compound (confirmed
        # by 'ditellanyl', same as `_diselenide.py`'s identical rule). A
        # bare (unprefixed) 'ditellanyl' name is the exception -- no
        # locant here either, so no parentheses (CID 173350364
        # 'ditellanylethane'), unlike the >=3-carbon case below where a
        # locant digit does sit next to it.
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


def _name_parent_chain(full_graph, carbon_graph, terminals, mol=None):
    chains = _longest_chains(carbon_graph)
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


def _name_benzene_ring_ditelluride_chain(mol, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since 'ditellanyl' has no suffix form (module
    docstring), a single, otherwise-unsubstituted benzene ring is senior
    to a chain of the same (plain-hydrocarbon) class regardless of the
    chain's length -- mirrors `_disulfide.py`'s
    `_name_benzene_ring_disulfide_chain`. The ring is always the parent
    hydride; the whole R-Te-Te- fragment hung *directly* off the ring is
    a single '...ditellanyl' substituent prefix on it (a chain spacer
    between the ring and the near tellurium is out of scope, same reason
    as `_disulfide.py`). Confirmed via PubChem: CID 101099279
    ('(methylditellanyl)benzene')."""
    te1_idx, te2_idx, c1, c2 = _validate_and_find_ditelluride(mol, aromatic_ring_atoms=ring_atoms)
    if c1 is None or c2 is None:
        raise UnsupportedStructure(
            "a -TeH terminal combined with a benzene-ring-substituent "
            "ditelluride is out of scope for this module"
        )
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter is not supported yet for a "
            "benzene-ring-substituent ditelluride"
        )

    full_graph = adjacency(mol)
    attachment = ring_chain_attachment(full_graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a ditelluride chain is not supported yet"
        )
    ring_atom, root = attachment

    if root == te1_idx:
        other_te, other_root = te2_idx, c2
    elif root == te2_idx:
        other_te, other_root = te1_idx, c1
    else:
        raise UnsupportedStructure(
            "a chain spacer between the benzene ring and the "
            "ditelluride's near tellurium is not supported yet (only a "
            "direct ring-tellurium bond is, see module docstring)"
        )

    sub_name, sub_compound = name_branch(full_graph, other_root, other_te, {}, mol=mol)
    if sub_compound:
        raise UnsupportedStructure(
            "a branched alkylditellanyl substituent is not supported yet"
        )
    return f"({sub_name}ditellanyl)benzene"


def name_ditelluride(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_ditelluride_chain(mol, ring_atoms)
    te1_idx, te2_idx, c1, c2 = _validate_and_find_ditelluride(mol)
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)

    if c1 is None or c2 is None:
        parent_root, parent_te = (c2, te2_idx) if c1 is None else (c1, te1_idx)
        terminals = {parent_te: _BARE_TERMINAL_NAME}
    else:
        size1 = len(_component_subgraph(carbon_graph, c1))
        size2 = len(_component_subgraph(carbon_graph, c2))

        if size1 == size2:
            _, compound_a = name_branch(full_graph, c1, te1_idx, {}, mol=mol)
            _, compound_b = name_branch(full_graph, c2, te2_idx, {}, mol=mol)
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

        sub_name, sub_compound = name_branch(full_graph, sub_root, sub_te, {}, mol=mol)
        if sub_compound:
            raise UnsupportedStructure(
                "a branched alkylditellanyl substituent is not supported yet"
            )
        terminals = {parent_te: sub_name + "ditellanyl"}

    parent_carbon_graph = _component_subgraph(carbon_graph, parent_root)
    chain, name = _name_parent_chain(full_graph, parent_carbon_graph, terminals, mol=mol)

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter must lie on the winning parent
    # chain; one on the ditellanyl (R') substituent branch is out of
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
