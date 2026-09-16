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

from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    component_subgraph,
    group_substituents,
    is_plain_benzene_ring,
    longest_chains,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_chain_forced_compound_terminals

_SULFUR = 16
_BARE_TERMINAL_NAME = "disulfanyl"


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


def _validate_and_find_disulfide(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom/ring/unsaturation rejections below so
    `name_disulfide`'s benzene-ring-substituent path (see
    `_name_benzene_ring_disulfide_chain`) can reuse this same validation
    for the rest of the molecule. Empty by default, so every other
    caller's behavior is unchanged."""
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
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if any(a not in aromatic_ring_atoms or b not in aromatic_ring_atoms for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see P-31 for "
            "alkenes/alkynes; not yet combined with a disulfide here)"
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

    s1, s2 = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == _SULFUR)
    c1 = None
    if s1.GetDegree() == 2:
        (c1,) = [n.GetIdx() for n in s1.GetNeighbors() if n.GetIdx() != s2.GetIdx()]
    c2 = None
    if s2.GetDegree() == 2:
        (c2,) = [n.GetIdx() for n in s2.GetNeighbors() if n.GetIdx() != s1.GetIdx()]
    return s1.GetIdx(), s2.GetIdx(), c1, c2


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
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
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
            substituents = substituents_for_chain_forced_compound_terminals(full_graph, candidate, terminals, mol=mol)
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_chain = key, name, candidate
    return best_chain, best_name


def _name_benzene_ring_disulfide_chain(mol, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since 'disulfanyl' has no suffix form (module
    docstring), a single, otherwise-unsubstituted benzene ring is senior
    to a chain of the same (plain-hydrocarbon) class regardless of the
    chain's length -- mirrors `_azide.py`'s
    `_name_benzene_ring_azide_chain`. The ring is always the parent
    hydride; the whole R-S-S- fragment hung directly off the ring is a
    single '...disulfanyl' substituent prefix on it. Confirmed via
    PubChem: CID 84234 ('(methyldisulfanyl)benzene'), CID 257711
    ('(ethyldisulfanyl)benzene').

    Deliberately narrower than `_azide.py`'s equivalent: only a *direct*
    ring-to-sulfur bond is supported here, not a chain spacer between the
    ring and the near sulfur (e.g. a benzyl disulfide, 'c1ccccc1CSSC') --
    PubChem's own name for that case ('(methyldisulfanyl)methylbenzene',
    CID 12779) nests the 'disulfanyl' prefix inside a further 'methyl'
    substituent layer in a way this project's general-purpose
    `name_branch` machinery (built for simpler terminal substituents like
    halogens) doesn't reproduce byte-for-byte, and resolving that
    discrepancy needs more research -- out of scope for this narrow first
    slice."""
    s1_idx, s2_idx, c1, c2 = _validate_and_find_disulfide(mol, aromatic_ring_atoms=ring_atoms)
    if c1 is None or c2 is None:
        raise UnsupportedStructure(
            "a perthiol (-S-SH) combined with a benzene-ring-substituent "
            "disulfide is out of scope for this module"
        )
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter is not supported yet for a "
            "benzene-ring-substituent disulfide"
        )

    full_graph = adjacency(mol)
    attachment = ring_chain_attachment(full_graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a disulfide chain is not supported yet"
        )
    ring_atom, root = attachment

    if root == s1_idx:
        other_sulfur, other_root = s2_idx, c2
    elif root == s2_idx:
        other_sulfur, other_root = s1_idx, c1
    else:
        raise UnsupportedStructure(
            "a chain spacer between the benzene ring and the disulfide's "
            "near sulfur is not supported yet (only a direct ring-sulfur "
            "bond is, see module docstring)"
        )

    sub_name, sub_compound = name_branch(full_graph, other_root, other_sulfur, {}, mol=mol)
    if sub_compound:
        raise UnsupportedStructure(
            "a branched alkyldisulfanyl substituent is not supported yet"
        )
    # Always parenthesized (P-16.3.3, same rule as the plain-chain path):
    # a bare, unprefixed 'disulfanyl' has no locant/prefix ambiguity to
    # guard against, but "[R']disulfanyl" always does.
    return f"({sub_name}disulfanyl)benzene"


def name_disulfide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_disulfide_chain(mol, ring_atoms)
    s1_idx, s2_idx, c1, c2 = _validate_and_find_disulfide(mol)
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)

    if c1 is None or c2 is None:
        parent_root, parent_s = (c2, s2_idx) if c1 is None else (c1, s1_idx)
        terminals = {parent_s: _BARE_TERMINAL_NAME}
    else:
        size1 = len(component_subgraph(carbon_graph, c1))
        size2 = len(component_subgraph(carbon_graph, c2))

        if size1 == size2:
            _, compound_a = name_branch(full_graph, c1, s1_idx, {}, mol=mol)
            _, compound_b = name_branch(full_graph, c2, s2_idx, {}, mol=mol)
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

        sub_name, sub_compound = name_branch(full_graph, sub_root, sub_s, {}, mol=mol)
        if sub_compound:
            raise UnsupportedStructure(
                "a branched alkyldisulfanyl substituent is not supported yet"
            )
        terminals = {parent_s: sub_name + "disulfanyl"}

    parent_carbon_graph = component_subgraph(carbon_graph, parent_root)
    chain, name = _name_parent_chain(full_graph, parent_carbon_graph, terminals, mol=mol)

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
