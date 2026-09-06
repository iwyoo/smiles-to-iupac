"""Naming of hydroperoxides (the '-peroxol' suffix, -OOH) on acyclic
saturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-56.1 (Chapter P-5, https://iupac.qmul.ac.uk/BlueBook/P5.html#5601):
  "It is now recommended to use the suffix 'peroxol' to provide
  substitutive names for hydroperoxides." Worked example: CH3-CH2-OOH ->
  'ethaneperoxol (PIN)' (replacing the older functional-class name 'ethyl
  hydroperoxide'). Structurally this is `_alcohol.py`'s '-ol' suffix with
  an extra bridging oxygen -- R-O-O-H rather than R-O-H.
- Like 'thiol' (`_thiol.py`), 'peroxol' begins with a consonant, so the
  parent hydride's final 'e' is never elided before it (P-16.3.3):
  'ethane' + 'peroxol' -> 'ethaneperoxol', not 'ethanperoxol'; by the same
  reasoning, a locant is inserted without eliding the 'e' either, e.g.
  'propane-1-peroxol' (mirroring `_thiol.py`'s 'propane-2-thiol', not
  `_alcohol.py`'s vowel-elided 'propan-2-ol'). The Blue Book's own text
  gives no worked example with a locant, but this follows directly from
  the same P-16.3.3 consonant rule already confirmed for 'thiol'.
- P-44.1.1 / P-44.4.1 / P-45.2 (Chapter P-4): same principal-chain and
  numbering machinery as `_alcohol.py`/`_thiol.py` -- the -OOH-bearing
  carbon gets the lowest available locant, ahead of substituent-prefix
  locants.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_alcohol.py`/`_thiol.py` apply here too (mononuclear parent, or a
  homogeneous two-carbon chain with exactly one substituent in total),
  e.g. 'ethaneperoxol' needs no locant.
- P-35.2.1 (Chapter P-3): halogen substituents are prefix-only and
  coexist freely with the 'peroxol' suffix, same as in `_alcohol.py`.

Scope, deliberately narrow (first pass at this functional group, mirroring
how `_thiol.py`/`_nitrile.py` each started in isolation): exactly one -OOH
group on an acyclic, saturated carbon chain, with no other heteroatom (in
particular no -OH, ether oxygen, or amine nitrogen) anywhere in the
molecule. Explicitly out of scope (raise `UnsupportedStructure`):
more than one -OOH (a bis-peroxol, P-56.1's multiplication -- future
work), any unsaturation (ene/yne coexistence), any ring, the P-56.2
chalcogen-analogue suffixes ('SO-thioperoxol' etc., a completely different
-S-OH/-Se-OH family), and any coexisting principal characteristic group
(Table 3.3 seniority competition, e.g. -OOH + -COOH) -- all separate,
unverified axes. One narrow exception to the "any ring" rule:
`_name_phenyl_chain_hydroperoxide` names a -OOH chain hanging off a
single plain, unsubstituted benzene ring (e.g.
'2-phenylethaneperoxol'), mirroring `_alcohol.py`/`_thiol.py`'s
identical benzene-ring-substituent path -- narrower than the acyclic
path: no chain unsaturation. The same function also covers the direct
case, -OOH with no chain at all bonded straight to the ring: unlike
'phenol' (`_alcohol.py`), P-56.1 has no retained ring-plus-suffix name
for this, so PubChem confirms benzene stays the parent with '-OOH' as a
plain 'hydroperoxy' prefix instead -- `OOc1ccccc1` -> "hydroperoxybenzene"
(same shape as `_nitro.py`'s 'nitrobenzene', not a suffix construction).
A substituted phenyl ring is still out of scope.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    non_single_bonds,
    ring_chain_attachment,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _hydroperoxide_oxygens(mol):
    """(attach, terminal) oxygen atoms of a plain -O-O-H hydroperoxide --
    two oxygens singly bonded to each other, one of degree 2 (bonded to a
    carbon) and the other of degree 1 with exactly one hydrogen -- or None
    if `mol` doesn't have exactly two oxygens shaped this way. Distinct
    from `_peroxide.py`'s R-O-O-R' shape, where both oxygens are degree 2."""
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(oxygens) != 2:
        return None
    o1, o2 = oxygens
    bond = mol.GetBondBetweenAtoms(o1.GetIdx(), o2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return None
    degrees = sorted((o1.GetDegree(), o2.GetDegree()))
    if degrees != [1, 2]:
        return None
    attach, terminal = (o1, o2) if o1.GetDegree() == 2 else (o2, o1)
    if terminal.GetTotalNumHs() != 1:
        return None
    (other,) = (n for n in attach.GetNeighbors() if n.GetIdx() != terminal.GetIdx())
    if other.GetAtomicNum() != 6:
        return None
    return attach, terminal


def has_hydroperoxide_shape(mol) -> bool:
    """True iff `mol` has exactly one plain -O-O-H hydroperoxide shape
    (P-56.1), regardless of whether the rest of the molecule is in scope."""
    return _hydroperoxide_oxygens(mol) is not None


def _validate_and_collect(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (site, exclude) where `site` is the carbon bearing -OOH and
    `exclude` is the {attach, terminal} oxygen indices to skip when
    collecting substituent branches off the parent chain.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_hydroperoxide`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_hydroperoxide`) can reuse
    this same validation for the rest of the molecule. Empty by default,
    so every other caller's behavior is unchanged."""
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a hydroperoxide's own two oxygens "
                "(P-56.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )

    oxygens = _hydroperoxide_oxygens(mol)
    if oxygens is None:
        raise UnsupportedStructure(
            "not a plain hydroperoxide (-OOH, P-56.1) shape; a peroxide "
            "R-O-O-R' is handled by a separate module (_peroxide.py)"
        )
    attach, terminal = oxygens
    (site,) = (n for n in attach.GetNeighbors() if n.GetIdx() != terminal.GetIdx())
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return site.GetIdx(), {attach.GetIdx(), terminal.GetIdx()}


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, locant, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "peroxol"

    if chain_length == 2 and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (here, the sole -OOH) in total omits the locant,
        # e.g. 'ethaneperoxol' (the Blue Book's own worked example).
        return alkane_name(2) + "peroxol"

    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)
    return f"{prefix}{stem}-{locant}-peroxol"


def _candidate_key(chain_length, locant, substituents):
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, locant, grouped)
    return (locant, -total_count, locant_set, citation_locants, name), name


def _substituents_for_chain(graph, chain, halogens, exclude):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in exclude]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_phenyl_chain_hydroperoxide(mol, ring_atoms):
    """Name a hydroperoxide whose -OOH lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 2-phenylethaneperoxol. The ring is
    cited as a 'phenyl' substituent prefix on the chain, which is the
    parent hydride, mirroring `_alcohol.py`'s `_name_phenyl_chain_alcohol`.
    Narrower than the acyclic path above: no chain unsaturation. Also
    handles the direct case (-OOH with no chain, bonded straight to the
    ring) as 'hydroperoxybenzene' -- see module docstring."""
    site, exclude = _validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "hydroperoxide chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain hydroperoxide is not supported yet"
        )
    ring_atom, chain_root = attachment
    if chain_root in exclude:
        # P-56.1 has no retained ring-plus-suffix name the way `_alcohol.py`'s
        # 'phenol' does for -OH; PubChem confirms benzene stays the parent
        # with '-OOH' cited as a plain 'hydroperoxy' prefix instead (same
        # shape as `_nitro.py`'s 'nitrobenzene', not a suffix construction).
        return "hydroperoxybenzene"
    chain, branches = longest_branched_chain_through(graph, site, ring_atoms, exclude)
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    halogens = halogen_substituents(mol)
    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        locant = position_of[site]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, locant, substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_acyclic_hydroperoxide(
    mol, site, exclude, extra_names=None, required_atoms=frozenset(), carbon_graph=None
):
    """`extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group demoted to a substituent prefix by
    `_seniority.senior_class` (e.g. a demoted amine's 'amino'), reused by
    `_hydroperoxide_amine.py` via `_coexisting_groups.py` so that module
    doesn't have to reimplement this chain search/candidate selection.
    `None` keeps the original halogens-only behavior unchanged.
    `required_atoms`: additional carbon atoms (e.g. a demoted amine's own
    carbon neighbor) that a candidate chain must also carry -- empty by
    default so existing callers are unaffected. `carbon_graph`: the
    carbon-only graph to search for the principal chain -- defaults to
    `carbon_adjacency(mol)` (unchanged behavior); `_ether_hydroperoxide.py`
    passes one with a coexisting ether's alkoxy-branch component already
    removed, the same reason and pattern as
    `_thiol.py`/`_ketone.py`/`_aldehyde.py`'s identical `carbon_graph`
    parameter (PR #428-430)."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **(extra_names or {})}
    chains = longest_chains(carbon_graph if carbon_graph is not None else carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = [chain for chain in chains if site in chain and required_atoms <= set(chain)]
    if not eligible:
        raise UnsupportedStructure(
            "the hydroperoxide-bearing carbon does not lie on a single "
            "longest carbon chain; a shorter principal chain capturing "
            "the -OOH group (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            locant = position_of[site]
            substituents = _substituents_for_chain(graph, candidate, halogens, exclude)
            key, name = _candidate_key(chain_length, locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_hydroperoxide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_hydroperoxide(mol, ring_atoms)
    site, exclude = _validate_and_collect(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see P-31 for "
            "alkenes/alkynes; not yet combined with a hydroperoxide here)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    return _name_acyclic_hydroperoxide(mol, site, exclude)
